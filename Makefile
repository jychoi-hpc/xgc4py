cybridge=xgc4py_c_bind
target=driver_xgc4py_c_bind

CC=gcc
CFLAGS= `python3-config --cflags`
LDFLAGS=`python3-config --ldflags` -L${HOME}/anaconda3/lib -lpython3.9 -Xlinker -rpath -Xlinker ${HOME}/anaconda3/lib

all:
	cython -3 $(cybridge).pyx
	$(CC) $(CFLAGS) -c $(target).c
	$(CC) $(CFLAGS) -c $(cybridge).c
	$(CC) *.o -o $(target) $(LDFLAGS)

clean:
	rm -f $(cybridge).{c,h,o} $(target).o $(target)
	rm -rf __pycache__

test:
	PYTHONHOME=$HOME/anaconda3 ./$(cybridge)
