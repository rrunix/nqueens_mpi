/**
MPI Version of the solution proposed by Jeff Somers (http://www.jsomers.com/nqueen_demo/nqueens.html)
*/

typedef unsigned long long SOLUTIONTYPE;
#include "nqmpi.h"

int board_size;

int depth;

SOLUTIONTYPE g_numsolutions;

int taskSent;

int taskReceived;

int size;

int rank;

void master_mpi_nqueen()
{

    int aQueenBitCol[MAX_BOARDSIZE];	/* marks colummns which already have queens */
    int aQueenBitPosDiag[MAX_BOARDSIZE];	/* marks "positive diagonals" which already have queens */
    int aQueenBitNegDiag[MAX_BOARDSIZE];	/* marks "negative diagonals" which already have queens */
    int aStack[MAX_BOARDSIZE + 2];	/* we use a stack instead of recursion */
    register int *pnStack;


    register int numrows = 0;	/* numrows redundant - could use stack */
    register unsigned int lsb;	/* least significant bit */
    register unsigned int bitfield;	/* bits which are set mark possible positions for a queen */
    int i;
    int odd = board_size & 1;	/* 0 if board_size even, 1 if odd */
    int board_minus = board_size - 1;	/* board size - 1 */
    int mask = (1 << board_size) - 1;	/* if board size is N, mask consists of N 1's */


    /* Initialize stack */
    aStack[0] = -1;		/* set sentinel -- signifies end of stack */

    /* NOTE: (board_size & 1) is true iff board_size is odd */
    /* We need to loop through 2x if board_size is odd */
    for (i = 0; i < (1 + odd); ++i) {
	/* We don't have to optimize this part; it ain't the 
	   critical loop */
	bitfield = 0;
	if (0 == i) {
	    /* Handle half of the board, except the middle
	       column. So if the board is 5 x 5, the first
	       row will be: 00011, since we're not worrying
	       about placing a queen in the center column (yet).
	     */
	    int half = board_size >> 1;	/* divide by two */
	    /* fill in rightmost 1's in bitfield for half of board_size
	       If board_size is 7, half of that is 3 (we're discarding the remainder)
	       and bitfield will be set to 111 in binary. */
	    bitfield = (1 << half) - 1;
	    pnStack = aStack + 1;	/* stack pointer */
	    aQueenBitCol[0] = aQueenBitPosDiag[0] = aQueenBitNegDiag[0] =
		0;
	} else {
	    /* Handle the middle column (of a odd-sized board).
	       Set middle column bit to 1, then set
	       half of next row.
	       So we're processing first row (one element) & half of next.
	       So if the board is 5 x 5, the first row will be: 00100, and
	       the next row will be 00011.
	     */
	    bitfield = 1 << (board_size >> 1);
	    numrows = 1;	/* prob. already 0 */

	    aQueenBitCol[0] = aQueenBitPosDiag[0] = aQueenBitNegDiag[0] =
		0;
	    aQueenBitCol[1] = bitfield;

	    /* Now do the next row.  Only set bits in half of it, because we'll
	       flip the results over the "Y-axis".  */
	    aQueenBitNegDiag[1] = (bitfield >> 1);
	    aQueenBitPosDiag[1] = (bitfield << 1);
	    pnStack = aStack + 1;	/* stack pointer */
	    *pnStack++ = 0;	/* we're done w/ this row -- only 1 element & we've done it */
	    bitfield = (bitfield - 1) >> 1;	/* bitfield -1 is all 1's to the left of the single 1 */
	}

	/* this is the critical loop */
	for (;;) {
	    /* could use 
	       lsb = bitfield ^ (bitfield & (bitfield -1)); 
	       to get first (least sig) "1" bit, but that's slower. */
	    lsb = -((signed) bitfield) & bitfield;	/* this assumes a 2's complement architecture */

	    // printf("printing, out: %d\n", numrows);
	    if (0 == bitfield) {
		bitfield = *--pnStack;	/* get prev. bitfield from stack */
		if (pnStack == aStack) {	/* if sentinel hit.... */
		    break;
		}
		--numrows;

		continue;
	    } else {
		bitfield &= ~lsb;	/* toggle off this bit so we don't try it again */
		if (numrows < board_minus) {	/* we still have more rows to process? */

		    if (numrows == depth) {
			int mpi_aQueenBitCol = aQueenBitCol[numrows] | lsb;

			int mpi_aQueenBitNegDiag =
			    (aQueenBitNegDiag[numrows] | lsb) >> 1;

			int mpi_aQueenBitPosDiag =
			    (aQueenBitPosDiag[numrows] | lsb) << 1;

			int mpi_bitfield = mask & ~(mpi_aQueenBitCol |
						    mpi_aQueenBitNegDiag |
						    mpi_aQueenBitPosDiag);

#ifdef NO_MPI
			g_numsolutions +=
			    slave_solver_mpi_nqueen(mpi_aQueenBitCol,
						    mpi_aQueenBitNegDiag,
						    mpi_aQueenBitPosDiag,
						    mpi_bitfield,
						    numrows + 1);
#else


			send_task(mpi_aQueenBitCol,
				  mpi_aQueenBitNegDiag,
				  mpi_aQueenBitPosDiag,
				  mpi_bitfield, numrows + 1);

#endif
			printf("\t M==> Called slave : sols now: %llu \n",
			       g_numsolutions);
			continue;
		    } else {
			int n = numrows++;
			aQueenBitCol[numrows] = aQueenBitCol[n] | lsb;
			aQueenBitNegDiag[numrows] =
			    (aQueenBitNegDiag[n] | lsb) >> 1;
			aQueenBitPosDiag[numrows] =
			    (aQueenBitPosDiag[n] | lsb) << 1;
			*pnStack++ = bitfield;
			/* We can't consider positions for the queen which are in the same
			   column, same positive diagonal, or same negative diagonal as another
			   queen already on the board. */
			bitfield =
			    mask & ~(aQueenBitCol[numrows] |
				     aQueenBitNegDiag[numrows] |
				     aQueenBitPosDiag[numrows]);
			continue;
		    }
		} else {
		    /* We have no more rows to process; we found a solution. */
		    /* Comment out the call to printtable in order to print the solutions as board position */
		    /* printtable(board_size, aQueenBitRes, g_numsolutions + 1);  */
		    ++g_numsolutions;
		    bitfield = *--pnStack;
		    --numrows;
		    continue;
		}
	    }
	}
    }
    waitUntilFinish();
    /* multiply solutions by two, to count mirror images */
    g_numsolutions *= 2;
}

SOLUTIONTYPE slave_solver_mpi_nqueen(int mpi_aQueenBitCol,
				     int mpi_aQueenBitNegDiag,
				     int mpi_aQueenBitPosDiag,
				     int mpi_bitfield, int mpi_numrows)
{
    SOLUTIONTYPE solutions = 0;

    int aQueenBitCol[MAX_BOARDSIZE];	/* marks colummns which already have queens */
    int aQueenBitPosDiag[MAX_BOARDSIZE];	/* marks "positive diagonals" which already have queens */
    int aQueenBitNegDiag[MAX_BOARDSIZE];	/* marks "negative diagonals" which already have queens */
    int aStack[MAX_BOARDSIZE + 2];	/* we use a stack instead of recursion */

    register int *pnStack;
    register int *pnStackStart;

    register int numrows = mpi_numrows;	/* numrows redundant - could use stack */
    register unsigned int lsb;	/* least significant bit */
    register unsigned int bitfield = mpi_bitfield;	/* bits which are set mark possible positions for a queen */
    int board_minus = board_size - 1;	/* board size - 1 */
    int mask = (1 << board_size) - 1;	/* if board size is N, mask consists of N 1's */

    /* Initialize stack */
    aStack[mpi_numrows - 1] = -1;	/* set sentinel -- signifies end of stack */
    aQueenBitCol[mpi_numrows] = mpi_aQueenBitCol;
    aQueenBitPosDiag[mpi_numrows] = mpi_aQueenBitPosDiag;
    aQueenBitNegDiag[mpi_numrows] = mpi_aQueenBitNegDiag;

    pnStack = aStack + mpi_numrows;
    pnStackStart = aStack + (mpi_numrows);
    *pnStack = bitfield;

    for (;;) {

	//lsb = bitfield ^ (bitfield & (bitfield -1)); slow way
	lsb = -((signed) bitfield) & bitfield;	/* this assumes a 2's complement architecture */
	if (bitfield == 0) {
	    if (pnStack == pnStackStart) {	/* if sentinel hit.... */
		break;
	    }
	    bitfield = *--pnStack;	/* get prev. bitfield from stack */
	    --numrows;
	    continue;
	}

	bitfield &= ~lsb;	/* toggle off this bit so we don't try it again */

	if (numrows < board_minus) {	/* we still have more rows to process? */
	    int n = numrows++;
	    aQueenBitCol[numrows] = aQueenBitCol[n] | lsb;
	    aQueenBitNegDiag[numrows] = (aQueenBitNegDiag[n] | lsb) >> 1;
	    aQueenBitPosDiag[numrows] = (aQueenBitPosDiag[n] | lsb) << 1;
	    *pnStack++ = bitfield;
	    /* We can't consider positions for the queen which are in the same
	       column, same positive diagonal, or same negative diagonal as another
	       queen already on the board. */
	    bitfield =
		mask & ~(aQueenBitCol[numrows] | aQueenBitNegDiag[numrows]
			 | aQueenBitPosDiag[numrows]);
	    continue;
	} else {
	    /* We have no more rows to process; we found a solution. */
	    ++solutions;
	    bitfield = *--pnStack;
	    --numrows;
	    continue;
	}
    }
    return solutions;
}

void waitUntilFinish()
{
    printf("Waiting slaves to finish: remaining %d from %d tasks\n",
	   taskReceived, taskSent);
    while (taskReceived < taskSent) {
	//Receive one
	MPI_Status info;
	unsigned long long sols;
	MPI_Recv(&sols, 1, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, TAG,
		 MPI_COMM_WORLD, &info);
	g_numsolutions += sols;
	printf("Slave %d finished\n", info.MPI_SOURCE);
	taskReceived++;
    }
    printf("All slaves finished\n");
}

void send_task(int mpi_aQueenBitCol, int mpi_aQueenBitNegDiag,
	       int mpi_aQueenBitPosDiag, int mpi_bitfield, int mpi_numrows)
{
    //Arguments to the slave
    int args[5] =
	{ mpi_aQueenBitCol, mpi_aQueenBitNegDiag, mpi_aQueenBitPosDiag,
	mpi_bitfield, mpi_numrows
    };
    if ((taskSent + 1) < size) {
	MPI_Send(&args[0], 5, MPI_INT, taskSent + 1, TAG, MPI_COMM_WORLD);
    } else {
	MPI_Request request;
	//Receive one
	MPI_Status info;
	unsigned long long sols;
	MPI_Recv(&sols, 1, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, TAG,
		 MPI_COMM_WORLD, &info);
	g_numsolutions += sols;
	//send other
	MPI_Isend(&args[0], 5,
		  MPI_INT, info.MPI_SOURCE, TAG, MPI_COMM_WORLD, &request);

	++taskReceived;
    }
    ++taskSent;
}

void handleMasterStart()
{
    time_t t1, t2;
    int i = 1;
    for (; i < size; i++) {
	MPI_Send(&board_size, 1, MPI_INT, i, TAG, MPI_COMM_WORLD);
    }
    time(&t1);
    printf("Start: \t %s\n\n\n", ctime(&t1));
    master_mpi_nqueen();
    time(&t2);
    printf("\nEnd: \t %s\n\n", ctime(&t2));
    tearing_down();
    printResults(&t1, &t2);
}

void handleSlaveStart()
{
    MPI_Status info;
    MPI_Recv(&board_size, 1, MPI_INT, MASTER_RANK, TAG, MPI_COMM_WORLD,
	     &info);
    while (1) {
	MPI_Request request;
	int args[5];
	MPI_Recv(&args[0], 5, MPI_INT, MASTER_RANK, TAG, MPI_COMM_WORLD,
		 &info);
	int length;
	MPI_Get_count(&info, MPI_INT, &length);
	if (length == 5) {
	    printf("\t S==> Slave start new task: %d %d %d %d %d \n",
		   args[0], args[1], args[2], args[3], args[4]);
	    SOLUTIONTYPE sols =
		slave_solver_mpi_nqueen(args[0], args[1], args[2], args[3],
					args[4]);
	    MPI_Isend(&sols, 1, MPI_LONG_LONG_INT, MASTER_RANK, TAG,
		      MPI_COMM_WORLD, &request);
	    printf("\t S==> Slave finish task\n");
	} else {
	    break;
	}
    }
}

void tearing_down()
{
    int i = 1;
    int none = 1;
    for (; i < size; i++) {
	MPI_Send(&none, 1, MPI_INT, i, TAG, MPI_COMM_WORLD);
    }
}

/* Print the results at the end of the run */
void printResults(time_t * pt1, time_t * pt2)
{

    printf
	("\n\n###########################################################################\n");
    double secs;
    int hours, mins, intsecs;

    if (g_numsolutions != 0) {
	printf("For board size %d, %llu solution%s found.\n", board_size,
	       g_numsolutions, (g_numsolutions == 1 ? "" : "s"));

    } else {
	printf("No solutions found.\n");
    }
    secs = difftime(*pt2, *pt1);
    intsecs = (int) secs;
    printf("Calculations took %d second%s.\n", intsecs,
	   (intsecs == 1 ? "" : "s"));

    /* Print hours, minutes, seconds */
    hours = intsecs / 3600;
    intsecs -= hours * 3600;
    mins = intsecs / 60;
    intsecs -= mins * 60;
    if (hours > 0 || mins > 0) {
	printf("Equals ");
	if (hours > 0) {
	    printf("%d hour%s, ", hours, (hours == 1) ? "" : "s");
	}
	if (mins > 0) {
	    printf("%d minute%s and ", mins, (mins == 1) ? "" : "s");
	}
	printf("%d second%s.\n", intsecs, (intsecs == 1 ? "" : "s"));
    }
    printf
	("###########################################################################\n\n");
}

/* main routine for N Queens program.*/
int main(int argc, char **argv)
{

    MPI_Status info;
#ifndef NO_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    if (rank == MASTER_RANK) {
	if (argc != 3) {
	    printf("Usage: nq <width of board> <server depth>\n");
	    return 0;
	}
	board_size = atoi(argv[1]);
	depth = atoi(argv[2]);

	handleMasterStart();

    } else {
	handleSlaveStart();
    }
    if (rank == MASTER_RANK) {
	printf("Tearing down master..\n");
    } else {
	printf("Tearing down slave %d..\n", rank);
    }
#ifndef NO_MPI
    MPI_Finalize();
#endif
    return 0;
}
