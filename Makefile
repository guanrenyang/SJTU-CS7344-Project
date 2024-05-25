# Makefile to compile and run MPI and serial ADI programs

# Compiler settings
MPICC=mpicc
CC=gcc
MPIRUN=mpirun
CFLAGS=-Wall -O2
MPICFLAGS=$(CFLAGS)
LDFLAGS=

# MPI program
MPI_SRC=1_adi.c
MPI_EXEC=adi

# Serial program
SERIAL_SRC=1_adi_serial.c
SERIAL_EXEC=adi_serial

# Build and run the MPI program
adi: $(MPI_SRC)
	$(MPICC) $(MPICFLAGS) -o $(MPI_EXEC) $(MPI_SRC) $(LDFLAGS)
	$(MPIRUN) -np 4 ./$(MPI_EXEC)

# Build and run the serial program
adi_serial: $(SERIAL_SRC)
	$(CC) $(CFLAGS) -o $(SERIAL_EXEC) $(SERIAL_SRC) $(LDFLAGS)
	./$(SERIAL_EXEC)

# Clean up build files
clean:
	rm -f $(MPI_EXEC) $(SERIAL_EXEC)
