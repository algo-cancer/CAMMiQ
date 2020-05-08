
if [ "$1" != "--cplex-dir" ]; then
	echo "Please clarify cplex directory!"
else
	cd src
	make download
	make divsufsort
	make CPLEXROOTDIR="../$2"
fi

