cd src
if [ -d "./parallel-divsufsort" ]; then
	echo "Destination path 'parallel-divsufsort' already exists."
else
	make downloads
fi
if [ -d "./robin-hood-hashing" ]; then
	echo "Destination path 'robin-hood-hashing' already exists."
else
	make downloadr
fi
DIVSUFSORTLIB=./parallel-divsufsort/lib/libdivsufsort.a
if [ -f "$DIVSUFSORTLIB" ]; then
	echo "Library 'libdivsufsort.a' already exists."
else
	make divsufsort
fi
if [ "$1" != "--cplex-dir" ]; then
	if [ "$1" != "--gurobi-dir" ]; then
		echo "Please clarify cplex or gurobi directory!"
	else
		if [ "$3" != "--gurobi-version" ]; then
			echo "Please clarify gurobi version!"
		else
			VERSION="$4"
			VERSION="${VERSION:0:3}"
			make GUROBIROOTDIR="$2" GRB_VERSION="${VERSION//.}" gurobi
		fi
	fi
else
	make CPLEXROOTDIR="$2" cplex
fi

