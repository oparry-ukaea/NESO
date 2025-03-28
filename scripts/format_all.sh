#!/usr/bin/env bash
REPO_ROOT=$( cd -- "$(realpath $( dirname -- "${BASH_SOURCE[0]}" )/..)" &> /dev/null && pwd )

clang_fmt_cmd="clang-format"
cmake_fmt_cmd="cmake-format"
python_fmt_cmd="black"

for cmd_ in $clang_fmt_cmd $cmake_fmt_cmd $python_fmt_cmd; do
    if ! command -v "$cmd_" &> /dev/null; then
        echo "$0: Command $cmd_ not found on your path"
        exit 2
    fi
done

if [[ -f "$REPO_ROOT/.clang-format" && -f "$REPO_ROOT/.cmake-format" ]]; then
    # clang-format
    find "$REPO_ROOT/src" "$REPO_ROOT/include" "$REPO_ROOT/test" "$REPO_ROOT/solvers" -iname \*.hpp -o -iname \*.cpp | xargs $clang_fmt_cmd -style=file:"$REPO_ROOT/.clang-format" -i
    # cmake-format
    cmake-format -c "$REPO_ROOT/.cmake-format" -i "$REPO_ROOT/CMakeLists.txt" # Do this one on it's own so we don't go into build etc
    find "$REPO_ROOT/src" "$REPO_ROOT/include" "$REPO_ROOT/test" "$REPO_ROOT/solvers" -iname CMakeLists.txt | xargs $cmake_fmt_cmd -c "$REPO_ROOT/.cmake-format" -i
    # black
    find "$REPO_ROOT/python" -iname \*.py | xargs $python_fmt_cmd
    exit 0;
else
    echo "ERROR: The files .clang-format and .cmake-format were not found in the repository root."
    exit 1;
fi

