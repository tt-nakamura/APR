NTL = -lntl -lgmp -L/usr/local/lib
APR = apr.o ZZlib.o Cycl.o

test: test.o $(APR)
	g++ test.o $(APR) $(NTL)