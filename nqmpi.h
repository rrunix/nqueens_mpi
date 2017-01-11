#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>


#define MIN_BOARDSIZE 2
#define MASTER_RANK 0
#define TAG 3
#define MAX_BOARDSIZE 26

//#define NO_MPI

typedef unsigned long long SOLUTIONTYPE;

SOLUTIONTYPE slave_solver_mpi_nqueen(int mpi_aQueenBitCol,
									 int mpi_aQueenBitNegDiag,
									 int mpi_aQueenBitPosDiag,
									 int mpi_bitfield, int mpi_numrows);

void send_task(int mpi_aQueenBitCol, int mpi_aQueenBitNegDiag,
			   int mpi_aQueenBitPosDiag, int mpi_bitfield,
			   int mpi_numrows);


void master_mpi_nqueen();


SOLUTIONTYPE slave_solver_mpi_nqueen(int mpi_aQueenBitCol,
									 int mpi_aQueenBitNegDiag,
									 int mpi_aQueenBitPosDiag,
									 int mpi_bitfield, int mpi_numrows);

void printResults(time_t* pt1, time_t* pt2);

void handleMasterStart();

void tearing_down();

void handleSlaveStart();

void waitUntilFinish();
