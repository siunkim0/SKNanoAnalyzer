# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SKNanoAnalyzer is a C++ physics analysis framework for CMS (Compact Muon Solenoid) experiment data analysis. It processes NanoAOD files for Run 2 (2016-2018) and Run 3 (2022-2023) data periods, providing tools for dilepton studies, ttbar reconstruction, and systematic uncertainty evaluations.

## Build Commands

```bash
# Initial setup (required once per session)
source setup.sh

# Build entire project
./scripts/build.sh

# Incremental rebuild
./scripts/rebuild.sh
```

## Running Analysis

```bash
# Submit analysis jobs
SKNano.py -a AnalyzerName -i 'SamplePattern*' -e 2022 -n 10

# Example: Run DiLepton analyzer on DYJets samples
SKNano.py -a DiLepton -i 'DYJets*' -e 2022 -n 50 --reduction 10

# Skimming mode for data reduction
SKNano.py -a Skim_AnalyzerName -i SampleName -e 2022 --skimming_mode
```

## Architecture Overview

### Core Components
- **DataFormats/**: Physics object classes (Electron, Muon, Jet, etc.) with analysis methods
- **AnalyzerTools/**: Analysis utilities, corrections, systematic helpers, ML integration  
- **Analyzers/**: Specific analysis implementations inheriting from AnalyzerCore
- **External dependencies**: KinematicFitter, RoccoR, jsonpog-integration, LHAPDF, libtorch

### Key Classes
- **AnalyzerCore**: Base class providing common functionality for all analyzers
- **SKNanoLoader**: Handles NanoAOD file loading and branch management
- **SystematicHelper**: Manages systematic uncertainty variations
- **MyCorrection**: Centralized correction and scale factor application

### Data Organization
- Era-based structure: `data/{era}/` contains corrections, efficiencies, and sample metadata
- Sample information in `data/{era}/Sample/CommonSampleInfo.json`
- Systematic configurations in YAML files (DiLeptonSystematic.yaml, VcbSystematic.yaml)

## Development Patterns

### Creating New Analyzers
1. Create header in `Analyzers/include/` inheriting from `AnalyzerCore`
2. Implement in `Analyzers/src/` with required methods:
   - `initializeAnalyzer()`
   - `executeEvent()`
   - Virtual destructor
3. Add to `Analyzers/include/AnalyzersLinkDef.hpp`
4. Use enum classes for constants instead of magic numbers/strings

### Systematic Uncertainties
- Use `SystematicHelper` class for systematic variations
- Implement variations using `variation` enum (nom, up, down)
- Apply corrections through `MyCorrection` class methods

### Physics Object Usage
```cpp
// Example pattern for lepton selection
RVec<Lepton *> leptons;
TString MuonID = "Tight"; // Example ID
for (auto &mu : *muons) {
    if (mu.Pass(MuonID)) {
        leptons.push_back(&mu);
    }
}
```

## Dependencies and Environment

### Required Packages (managed via conda/mamba)
- ROOT 6.32+
- correctionlib (CMS corrections)
- onnxruntime-cpp (ML inference)
- PyTorch/libtorch (deep learning)
- LHAPDF (auto-installed via scripts/install\_lhapdf.sh)

### Configuration
- Personal config: `config/config.$USER` (system type, package manager, notifications)
- Sample lists: `SampleLists/{era}/` for different analysis categories

### Submodules (CERN gitlab)

Two submodules are hosted on `gitlab.cern.ch` and require **CERN gitlab SSH access** to clone properly:
- `external/RoccoR` — Rochester muon momentum corrections (`ssh://git@gitlab.cern.ch:7999/hyseo/RoccoR.git`)
- `external/jsonpog-integration` — CMS POG correction JSON files (`ssh://git@gitlab.cern.ch:7999/cms-nanoAOD/jsonpog-integration.git`)

**With CERN gitlab access** (normal setup):
```bash
git submodule update --init external/RoccoR external/jsonpog-integration
```

**Without CERN gitlab access** (workaround so `git status` / `git commit` don't fail with `fatal: Not a git repository: external/RoccoR/...`): create empty stub gitdirs so git's submodule machinery is satisfied. The stubs live only in `.git/` and are not pushed.

```bash
# Create stub for RoccoR
mkdir -p .git/modules/external/RoccoR/{objects/info,objects/pack,refs/heads,refs/tags}
printf "ref: refs/heads/master\n" > .git/modules/external/RoccoR/HEAD
printf "[core]\n\trepositoryformatversion = 0\n\tfilemode = true\n\tbare = false\n\tworktree = ../../../../external/RoccoR\n" > .git/modules/external/RoccoR/config

# Create stub for jsonpog-integration (note: this gitdir lives under .git/worktrees/... due to historical layout)
mkdir -p .git/worktrees/SKNanoAnalyzer-v13/modules/external/jsonpog-integration/{objects/info,objects/pack,refs/heads,refs/tags}
printf "ref: refs/heads/master\n" > .git/worktrees/SKNanoAnalyzer-v13/modules/external/jsonpog-integration/HEAD
printf "[core]\n\trepositoryformatversion = 0\n\tfilemode = true\n\tbare = false\n\tworktree = ../../../../../../external/jsonpog-integration\n" > .git/worktrees/SKNanoAnalyzer-v13/modules/external/jsonpog-integration/config
```

**Switching from stubs to real submodules** (once CERN access is granted):
```bash
rm -rf .git/modules/external/RoccoR .git/worktrees/SKNanoAnalyzer-v13
git submodule update --init external/RoccoR external/jsonpog-integration
```

The stub workaround unblocks committing to the main repo, but the submodule contents themselves remain unusable (no corrections data) until proper init.

## Job Management

### HTCondor Integration
- Jobs submitted via HTCondor with DAG workflows
- Configurable memory (default 2GB) and CPU resources
- Automatic output file management and merging
- Telegram bot notifications (optional)

### Sample Selection
- Regex patterns for flexible sample matching
- Era-specific sample organization (2016preVFP, 2017, 2018, 2022, 2022EE, 2023, etc.)
- Cross-section and luminosity handled automatically

## Multi-Era Support

The framework supports unified analysis across data-taking periods:
- **Run 2**: 2016preVFP, 2016postVFP, 2017, 2018
- **Run 3**: 2022, 2022EE, 2023, 2023BPix

Era-specific corrections and configurations are automatically loaded based on the `-e` flag in job submission.
