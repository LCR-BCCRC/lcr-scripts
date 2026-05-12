#!/bin/bash
set -euo pipefail

python -c "import SigProfilerExtractor; import SigProfilerMatrixGenerator; import sigProfilerPlotting"
