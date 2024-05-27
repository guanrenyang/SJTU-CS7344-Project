# Makefile to compile and run MPI and serial ADI programs as well as Floyd's algorithm programs

# Compiler settings
MPICC=mpicc
CC=gcc
MPIRUN=mpirun
CFLAGS=-Wall -O2
MPICFLAGS=$(CFLAGS)
LDFLAGS=

# MPI program for ADI
MPI_SRC=1_adi.c
MPI_EXEC=adi

# Serial program for ADI
SERIAL_SRC=1_adi_serial.c
SERIAL_EXEC=adi_serial

# MPI program for Floyd's algorithm
FLOYD_MPI_SRC=3_floyd.c
FLOYD_MPI_EXEC=floyd

# Serial program for Floyd's algorithm
FLOYD_SERIAL_SRC=3_floyd_serial.c
FLOYD_SERIAL_EXEC=floyd_serial

# Build and run the MPI program for ADI
adi: $(MPI_SRC)
	$(MPICC) $(MPICFLAGS) -o $(MPI_EXEC) $(MPI_SRC) $(LDFLAGS)
	$(MPIRUN) -np 4 ./$(MPI_EXEC)

# Build and run the serial program for ADI
adi_serial: $(SERIAL_SRC)
	$(CC) $(CFLAGS) -o $(SERIAL_EXEC) $(SERIAL_SRC) $(LDFLAGS)
	./$(SERIAL_EXEC)

# Build and run the MPI program for Floyd's algorithm
floyd: $(FLOYD_MPI_SRC)
	$(MPICC) $(MPICFLAGS) -o $(FLOYD_MPI_EXEC) $(FLOYD_MPI_SRC) $(LDFLAGS)
	$(MPIRUN) -np 4 ./$(FLOYD_MPI_EXEC)

# Build and run the serial program for Floyd's algorithm
floyd_serial: $(FLOYD_SERIAL_SRC)
	$(CC) $(CFLAGS) -o $(FLOYD_SERIAL_EXEC) $(FLOYD_SERIAL_SRC) $(LDFLAGS)
	./$(FLOYD_SERIAL_EXEC)

# Clean up build files
clean:
	rm -f $(MPI_EXEC) $(SERIAL_EXEC) $(FLOYD_MPI_EXEC) $(FLOYD_SERIAL_EXEC)
