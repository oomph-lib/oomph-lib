setup:
	echo "Can't use this yet. WIP. Need to make sure it runs safely."
	mkdir build && cd build
	if command -v ninja &> /dev/null; then
		cmake -G Ninja ..
		ninja && ninja install
	else
		cmake ..
		make && make install
	fi

check:
	cd demo_drivers
	mkdir build
	cd build
	if command -v ninja &> /dev/null; then
		cmake -G Ninja ..
	else
		cmake ..
	fi
	ctest -j4

cmake-format:
	set -f
	find . \( -name '*.cmake' -o -name 'CMakeLists.txt' \) -exec cmake-format -c .cmake-format.json -i {} \;
	set +f

clang-format:
	echo "Now starting clang-format run on src/ and demo_drivers/ folder..."
	./scripts/run_clang_format_on_all.sh ./src
	./scripts/run_clang_format_on_all.sh ./demo_drivers
	echo "Finished running clang-format."