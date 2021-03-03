parallel “echo {} && gunzip -c {} | wc -l | awk ‘{d=\$1; print d/4;}” ::: *.gz
