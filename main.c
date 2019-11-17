#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void usr_input(double ***A, unsigned* psize, double** b);

void res_output(double** A,unsigned size,double* b,double* x);

void gauss_forward(double** A,unsigned size,double* b);

void gauss_back(double** A,unsigned size,double* b,double* x);

int main (int argc,char* argv[])
 {
 double **A,*b;
 unsigned size;

 usr_input(&A,&size,&b);

 double* x=(double*) malloc(sizeof(double)*size);

 gauss_forward(A,size,b);
 gauss_back(A,size,b,x);

// memset(x,0.,size);

 res_output(A,size,b,x);
 for(unsigned i=0;i<size;++i)
  free(A[i]);
 free(A);
 free(b);
 free(x);
 }

void usr_input(double*** A, unsigned *psize, double** b)
 {
 printf("Input size of matrix A: ");
 scanf("%i",psize);

 *A=(double**) malloc(sizeof(double*)*(*psize));
 for(unsigned i=0;i<*psize;++i)
  (*A)[i]=(double*) malloc(sizeof(double)*(*psize));
 *b=(double*) malloc(sizeof(double)*(*psize));

 for(unsigned i=0;i<*psize;++i)
  {
  printf("Fill matrix string A[%i]: ",i);
  fflush(stdout);
  for(unsigned j=0;j<*psize;++j)
   scanf("%lf",(*A)[i]+j);
  }

// printf("*(A[i]+j)=%0.3lf\n",*((*A)[1]+2));
// printf("A[1][2]=%0.3lf\n",(*A)[1][2]);

 printf("Fill vector b: ");
  for(unsigned j=0;j<*psize;++j)
   scanf("%lf",*b+j);
 }

void res_output(double** A,unsigned size,double* b,double* x)
 {
 printf("\nThe solution of linear system:\n");
 for(unsigned i=0;i<size;++i)
  {
  for(unsigned j=0;j<size;++j)
   printf("%.3lf\t",A[i][j]);
  printf("| %.3lf\n",b[i]);
  }

 printf("is:\n");
 for(unsigned i=0;i<size;++i)
  printf("%.3lf\n",x[i]);
 }

void gauss_forward(double** A,unsigned size,double* b)
 {
 unsigned idx_max;
 double *B,coef,c;
//Цикл #1 Проход по всем строками
 for(unsigned i=0;i<size;++i)
  {
//Поиск максимума в столбце
  idx_max=i;
//Цикл #1.1 Проход по нижним строкам (Можно распараллелить по данным и в конце MPI_Reduce(,MPI_MAX,);)
  for(unsigned k=i+1;k<size;++k)
   if(fabs(A[idx_max][i])<fabs(A[k][i]))
    idx_max=k;

//Перестановка строк (Быстрей будет всем процессам переставить строки (недолго))
  if(idx_max!=i)
   {
   B=A[i];
   A[i]=A[idx_max];
   A[idx_max]=B;
   c=b[i];
   b[i]=b[idx_max];
   b[idx_max]=c;
   }

//Ранг меньше размерности (главный процесс может вывести сообщение)
  if(A[i][i]==0)
   {
   printf("Rank of A < size(=%i)\n",size);
   exit(1);
   }
//Цикл #1.2 Проход по нижним строкам (можно распараллелить по строкам ЛИБО ЗДЕСЬ
//(надо подумать над равномерностью - норм, если число процессов | число итераций.
//Собираем, обмениваясь строками(не очень),наверно, через MPI_Allgather();)
  for(unsigned k=i+1;k<size;++k)
   {
//  if(A[k*size+i]!=0) оптимизация, потом можно проверить
   coef=(A[k][i]/A[i][i]);
//Цикл #1.2.1 Проход по правым координатам (ЛИБО ЗДЕСЬ, проблема, когда
//число процессов !| число итераций, выражена так же. Собираем,обмениваясь координатами,
//наверно, через MPI_Allgather();
//(скорее всего так же, а то и хуже, т.к. число пересылок больше в o(n) раз))
   for(unsigned j=i;j<size;++j)
    A[k][j]-=coef*A[i][j];
   b[k]-=coef*b[i];
   }
  }
 }

void gauss_back(double** A,unsigned size,double* b,double* x)
 {
//Цикл #2 Проход по всем строкам
 for(unsigned i=size-1;i<size;--i) //т.к. тип беззнаковый
  {
  x[i]=b[i];
//Цикл #2.1 Проход по правым координатам (можно распараллелить по координатам здесь.
//Собираем MPI_Reduce(,MPI_SUM,);)
  for(unsigned j=i+1;j<size;++j)
   x[i]-=A[i][j]*x[j];
  x[i]/=A[i][i];
  }
 }
