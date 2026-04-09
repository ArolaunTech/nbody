sh build.sh
if [ $? -eq 0 ]; then
    ./nbody
else
	echo "Build failed"
fi