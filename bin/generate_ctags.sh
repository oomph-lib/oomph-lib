
# Directories to include
dir="demo_drivers/ src/ self_test/ user_drivers/ user_src/"

# Find the c++ files in those directories (and not in .ccls-cache)
find $dir -type f \( -name "*.h" -o -name "*.hpp" -o -name "*.cc" -o -name "*.cpp" \) -not -path "*/.ccls-cache/*" > cpp_files.txt

# Tag those files
ctags --c++-kinds=+p --fields=+iaS --extras=+q --extras=f --languages=C++ -L cpp_files.txt --verbose
