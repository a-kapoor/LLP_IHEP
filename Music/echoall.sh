echo "---------------------------------"
echo "Python Version"
echo $(python -c "import sys; print(sys.version_info);" || python -c "import sys; print sys.version_info;")
echo $(python -c "import sys; print(sys.version);" || python -c "import sys; print sys.version;")
echo "---------------------------------"
