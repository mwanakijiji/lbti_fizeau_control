# Fizeau control scripts

To rsync the .py files only, run in this directory (as of 2019 Dec 27)
rsync -vr --include="*/" --include="*.py" --exclude="*" observer@lbti-data:~/home/observer/scripts_es/lbti_fizeau_control/* .
