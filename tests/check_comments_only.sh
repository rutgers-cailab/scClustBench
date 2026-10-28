#!/bin/bash
# Confirms that edits between two copies of the repository changed comments only:
# R files must parse to identical expressions, Python files to identical ASTs,
# and shell scripts must match once full-line comments and blank lines are removed.
# Usage: bash tests/check_comments_only.sh <old copy> [<new copy, default .>]
set -uo pipefail
OLD=$1; NEW=${2:-.}; fail=0
while IFS= read -r f; do
  [ -f "$NEW/$f" ] || { echo "MISSING $f"; fail=1; continue; }
  case $f in
    *.R) Rscript -e 'a <- commandArgs(TRUE); d <- function(p) deparse(parse(p, keep.source = FALSE));
                     if (!identical(d(a[1]), d(a[2]))) quit(status = 1)' "$OLD/$f" "$NEW/$f" ;;
    *.py) python3 -c 'import ast,sys; d=lambda p: ast.dump(ast.parse(open(p).read())); sys.exit(d(sys.argv[1]) != d(sys.argv[2]))' "$OLD/$f" "$NEW/$f" ;;
    *.sh|*.sbatch) cmp -s <(grep -vE '^\s*(#[^!]|#$|$)' "$OLD/$f") <(grep -vE '^\s*(#[^!]|#$|$)' "$NEW/$f") ;;
    *) true ;;
  esac || { echo "CODE CHANGED $f"; fail=1; }
done < <(cd "$OLD" && find . -type f \( -name '*.R' -o -name '*.py' -o -name '*.sh' -o -name '*.sbatch' \) | sed 's|^\./||')
[ $fail = 0 ] && echo "comments-only check passed"
exit $fail
