// Compare two SKNanoAnalyzer output files (histograms + trees) bin-by-bin and
// entry-by-entry. Prints IDENTICAL if everything matches, DIFFER otherwise.
//
// Usage:
//   rot && root -l -b -q 'scripts/compare_outputs.C("before.root","after.root")'

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH1.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TString.h"

static long g_nCompared = 0;
static long g_nDiff = 0;

static void report(const TString &what, const TString &detail) {
    if (g_nDiff < 50) std::cout << "  DIFF: " << what << " : " << detail << std::endl;
    g_nDiff++;
}

static void compareHist(TH1 *h1, TH1 *h2, const TString &path) {
    g_nCompared++;
    if (h1->GetNcells() != h2->GetNcells()) {
        report(path, Form("ncells %d vs %d", h1->GetNcells(), h2->GetNcells()));
        return;
    }
    if (h1->GetEntries() != h2->GetEntries())
        report(path, Form("entries %g vs %g", h1->GetEntries(), h2->GetEntries()));
    for (int i = 0; i < h1->GetNcells(); i++) {
        if (h1->GetBinContent(i) != h2->GetBinContent(i)) {
            report(path, Form("bin %d content %.17g vs %.17g", i, h1->GetBinContent(i), h2->GetBinContent(i)));
            break;
        }
        if (h1->GetBinError(i) != h2->GetBinError(i)) {
            report(path, Form("bin %d error %.17g vs %.17g", i, h1->GetBinError(i), h2->GetBinError(i)));
            break;
        }
    }
}

static void compareTree(TTree *t1, TTree *t2, const TString &path) {
    g_nCompared++;
    const Long64_t n1 = t1->GetEntries();
    const Long64_t n2 = t2->GetEntries();
    if (n1 != n2) {
        report(path, Form("tree entries %lld vs %lld", n1, n2));
        return;
    }
    std::vector<TString> leafNames;
    for (auto obj : *(t1->GetListOfLeaves())) leafNames.push_back(obj->GetName());
    if ((size_t)t2->GetListOfLeaves()->GetEntries() != leafNames.size())
        report(path, "different number of leaves");
    for (Long64_t i = 0; i < n1; i++) {
        t1->GetEntry(i);
        t2->GetEntry(i);
        for (const auto &lname : leafNames) {
            TLeaf *l1 = t1->GetLeaf(lname);
            TLeaf *l2 = t2->GetLeaf(lname);
            if (!l1 || !l2) { report(path + "/" + lname, "leaf missing"); continue; }
            const int len1 = l1->GetLen(), len2 = l2->GetLen();
            if (len1 != len2) {
                report(path + "/" + lname, Form("entry %lld len %d vs %d", i, len1, len2));
                continue;
            }
            for (int k = 0; k < len1; k++) {
                const double v1 = l1->GetValue(k), v2 = l2->GetValue(k);
                if (v1 != v2 && !(std::isnan(v1) && std::isnan(v2))) {
                    report(path + "/" + lname, Form("entry %lld[%d] %.17g vs %.17g", i, k, v1, v2));
                    break;
                }
            }
        }
        if (g_nDiff > 1000) { std::cout << "  ... too many diffs, aborting tree " << path << std::endl; return; }
    }
}

static void compareDir(TDirectory *d1, TDirectory *d2, const TString &prefix) {
    TIter next(d1->GetListOfKeys());
    while (TKey *key = (TKey *)next()) {
        const TString name = key->GetName();
        const TString path = prefix.Length() ? prefix + "/" + name : name;
        TObject *o1 = d1->Get(name);
        TObject *o2 = d2->Get(name);
        if (!o2) { report(path, "missing in second file"); continue; }
        if (o1->InheritsFrom(TDirectory::Class())) {
            compareDir((TDirectory *)o1, (TDirectory *)o2, path);
        } else if (o1->InheritsFrom(TTree::Class())) {
            compareTree((TTree *)o1, (TTree *)o2, path);
        } else if (o1->InheritsFrom(TH1::Class())) {
            compareHist((TH1 *)o1, (TH1 *)o2, path);
        }
    }
    // objects present only in the second file
    TIter next2(d2->GetListOfKeys());
    while (TKey *key = (TKey *)next2()) {
        if (!d1->Get(key->GetName()))
            report(prefix.Length() ? prefix + "/" + key->GetName() : TString(key->GetName()),
                   "missing in first file");
    }
}

void compare_outputs(const char *file1, const char *file2) {
    TFile *f1 = TFile::Open(file1);
    TFile *f2 = TFile::Open(file2);
    if (!f1 || f1->IsZombie() || !f2 || f2->IsZombie()) {
        std::cout << "ERROR: cannot open input files" << std::endl;
        return;
    }
    compareDir(f1, f2, "");
    std::cout << "Compared " << g_nCompared << " objects, " << g_nDiff << " difference(s)" << std::endl;
    std::cout << (g_nDiff == 0 ? "IDENTICAL" : "DIFFER") << std::endl;
}
