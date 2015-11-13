/* ============================================================================

   Copyright (C) 1997, 2007, 2010  Konrad Bernloehr 

   This file is part of the eventio/hessio library.

   The eventio/hessio library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this library. If not, see <http://www.gnu.org/licenses/>.

============================================================================ */

/** @file dhsort.c
 *  @short dhsort - double type number heapsort 
 *
 *  @author  Konrad Bernloehr 
 *  @date    $Date: 2010/07/20 13:37:45 $
 *  @version $Revision: 1.6 $

\verbatim

            Based on algorithms by Jon Bentley [Communications of the ACM v
            28 n 3 p 245 (Mar 85) and v 28 n 5 p 456 (May 85)], and the
            sort interface routines by Allen I.  Holub [Dr.  Dobb's Journal
            #102 (Apr 85)].

    Notes...
            This routine sorts N doubles in worst-case time proportional to
            N*log(N).  The heapsort was discovered by J.  W.  J.  Williams
            [Communications of the ACM v 7 p 347-348 (1964)] and is
            discussed by D.  E.  Knuth [The Art of Computer Programming,
            Volume 3: Sorting and Searching, Addison-Wesley, Reading,
            Mass., 1973, section 5.2.3].

            This algorithm depends on a portion of an array having the
            "heap" property.  The array X has the property heap[L,U] if:

        	    for all	 L, i, and U
        	    such that	 2L <= i <= U
        	    we have	 X[i div 2] <= X[i]
\endverbatim
*/
/* ================================================================ */

#include "initial.h"

#include "dhsort.h"

/** Perform a heap sort on a double array starting at dnum. */

void dhsort(double *dnum, int nel)
{
   int i=0, c=0, L=0;
   double *DNum;
   double DTmp;

   DNum=dnum-1;
   for (i=nel/2; i>=1; i--)
   {
      /* Enlarge heap */
      L = i;
      while(1)
      {
         c = L+L;
	 if (c>nel) break;
	 if ( (c < nel) && (DNum[c+1]>DNum[c]) ) c++;
	 if (DNum[L] >= DNum[c]) break;
#ifdef DEBUG_DHSORT
 printf("Exchange %f and %f (pos %d and %d)\n",DNum[L],DNum[c],L,c);
#endif
	 DTmp = DNum[L];
	 DNum[L] = DNum[c];
	 DNum[c] = DTmp;
	 L = c;
      }
   }
   
   for (i=nel; i>=2; )
   {
#ifdef DEBUG_DHSORT
 printf("Exchange %f and %f (pos 1 and %d)\n",*dnum,DNum[i],i);
#endif
      DTmp = *dnum;
      *dnum = DNum[i];
      DNum[i] = DTmp;
      i--;
      /* Enlarge heap */
      L = 1;
      while(1)
      {
      	 c = L+L;
	 if (c>i) break;
	 if ( (c < i) && (DNum[c+1]>DNum[c]) ) c++;
	 if (DNum[L] >= DNum[c]) break;
#ifdef DEBUG_DHSORT
 printf("Exchange %f and %f (pos %d and %d)\n",DNum[L],DNum[c],L,c);
#endif
	 DTmp = DNum[L];
	 DNum[L] = DNum[c];
	 DNum[c] = DTmp;
	 L = c;
      }
   }

}


#ifdef DEBUG_DHSORT

int main (int argc, char **argv)
{
   double d[1000], dt;
   int i, nd = 0;
   
   for (i=1; i<argc && i<1000; i++)
   {
      if ( sscanf(argv[i],"%lf",&dt) == 1 )
         d[nd++] = dt;
   }
   printf("Before dhsort:\n");
   for (i=0; i<nd; i++)
       printf("%f\n",d[i]);
   dhsort(d,nd);
   printf("After dhsort:\n");
   for (i=0; i<nd; i++)
       printf("%f\n",d[i]);
   return 0;
}

#endif
