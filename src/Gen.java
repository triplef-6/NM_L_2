import static java.lang.Math.abs;

/**
 * Hello world!
 *
 */
public class Gen
{
    public Gen() {
    }

    /**
     * показать матрицу
     * @param a матрица
     * @param n размер матрицы
     */
    public void print_matr(double[][] a, int n) {
         System.out.println(" ");
         for (int i = 0 ; i < n; i++) {
            for (int j = 0; j < n; j++) {
                System.out.printf("%5.2f ", a[i][j]);
            }
            System.out.println();
        }
    }

    /**
     * кубическая норма матрицы
     * @param a матрица
     * @param n размер матрицы
     * @return норма
     */
    public double matr_inf_norm(double[][] a, int n) {
        double norm = 0.;

        for (int i = 0; i < n; i++) {
            double s = 0.;
            for (int j = 0; j < n; j++) {
                s += Math.abs( a[i][j] );
            }
            if (s > norm) {
                norm = s;
            }
        }

        return norm;
    }

    /**
     * умножение матриц
     * @param a матрица 1
     * @param b матрица 2
     * @param n размер матрицы
     * @return произведение матриц
     */
    public double[][] matr_mul(double[][] a,  double[][] b, int n) {
        double[][] c = new double[n][n];
        for(int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int k = 0;
                for (c[i][j] = 0.; k < n; k++) {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return c;
    }

    /**
     * ?поиск ортогональной матрицы
     * @param Q
     * @param n
     */
    public void Q_matrix(double[][] Q, int n) {
        double  q;

        double curr, next = 1.;
        for (int j=0; j<n-1; j++) {
            curr = next;
            next += 1.;
            q = 1. / Math.sqrt(curr * next);
            for (int i = 0; i <= j; i++) {
                Q[i][j] = q;
            }
            Q[j+1][j] = -Math.sqrt(curr / next);
            for(int i = j + 2; i < n; i++) {
                Q[i][j] = 0.;
            }
        }
        q = 1. / Math.sqrt(n);
        for (int i = 0; i < n; i++) {
            Q[i][n - 1] = q;
        }
    }

    /**
     * генератор матриц
     * @param a генерируемая матрица
     * @param a_inv обратная матрица
     * @param n размерность
     * @param alpha альфа
     * @param beta бетта
     * @param sign_law параметры моделирования
     * @param lambda_law параметры моделирования
     * @param variant параметры моделирования
     * @param schema параметры моделирования
     * 1 2 0 1 симметричная
     * 1 2 1 1 простой структуры
     * 0 0 2 1 жордановы клетки
     */
    public GenMatrixParam mygen(double[][] a, double[][] a_inv, int n, double alpha, double beta, int sign_law, int lambda_law, int variant, int schema, boolean print_flag) {
        int i,j,k;
        GenMatrixParam genMatrixParam = new GenMatrixParam();

//        System.out.println("   M A T R I X  G E N.  ");
//        System.out.println("              N = "+ n);
//        System.out.println(" | lambda_min | = "+alpha);
//        System.out.println(" | lambda_max | = "+beta);

        double[] lambda = new double[n];

        // распределение знаков
//        System.out.println(" sign_law = "+sign_law);

        double[] sign   = new double[n];
        for( i=0; i<n; i++ ) sign[i] = 1.;

        switch(sign_law)
        {
            case -1:
                for( i=0; i<n; i++ ) sign[i] =-1.;
                break;

            case 0:
                sign[0] = 1.;
                for( i=1; i<n; i++ ) sign[i] =-sign[i-1];
                break;

            //другие законы распределения знаков
            // ...

        }
/*	for( i=0; i<n; i++ ) cout<<sign[i]<<" ";
	cout<<endl;
*/

        //распределение собственнных чисел
//        System.out.println(" lambda_law = "+lambda_law);

        double[] kappa   = new double[n];
        for( i=0; i<n; i++ ) kappa[i] = (double)i/(double)(n-1);
        switch(lambda_law)
        {
            case 1:
                System.out.println(" kappa = sqrt( ) ");
                for( i=0; i<n; i++ ) kappa[i] = Math.sqrt( kappa[i] );
                break;

            case 2:
//                System.out.println(" kappa = sin( ) ");
                double pi_half = Math.acos(-1.)*0.5;
                for( i=0; i<n; i++ ) kappa[i] = Math.sin( pi_half*kappa[i] );
                break;

            //другие законы распределения собственных чисел
            // ...

        }
/*	for( i=0; i<n; i++ ) cout<<kappa[i]<<" ";
	cout<<endl;
*/


        double[] J = new double[n];
        for( i=0; i<n; i++ ) J[i] = sign[i]*( ( 1.-kappa[i] )*alpha + kappa[i]*beta );
/*	for( i=0; i<n; i++ ) cout<<J[i]<<" ";
	cout<<endl;
*/

        double[] J_inv = new double[n];
        for( i=0; i<n; i++ ) J_inv[i] = 1./J[i];

        double[][] Q = new double [n][];
        for( i=0; i<n; i++ ) Q[i] = new double [n];

        double[] aa = new double[3];


//        System.out.println(" variant = "+variant);

        switch( variant )
        {
            case 0: //симметричная матрица
//                System.out.println(" simmetric matrix:");
//                System.out.println(" schema = "+schema);
                switch( schema )
                {
                    case 1:
                        Q_matrix ( Q, n);

                        for( a[0][0]=0., k=0; k<n; k++ ) a[0][0] += Q[0][k]*J[k]*Q[0][k];
                        for( j=1; j<n; j++ )
                        {
                            for( a[0][j]=0., k=j-1; k<n; k++ ) a[0][j] += Q[0][k]*J[k]*Q[j][k];
                            a[j][0] = a[0][j];
                        }
                        for( i=1; i<n; i++ )
                        {
                            for( a[i][i]=0., k=i-1; k<n; k++ ) a[i][i] += Q[i][k]*J[k]*Q[i][k];
                            for( j=i+1; j<n; j++ )
                            {
                                for( a[i][j]=0., k=j-1; k<n; k++ ) a[i][j] += Q[i][k]*J[k]*Q[j][k];
                                a[j][i] = a[i][j];
                            }
                        }

//_______
                        for( a_inv[0][0]=0., k=0; k<n; k++ ) a_inv[0][0] += Q[0][k]*J_inv[k]*Q[0][k];
                        for( j=1; j<n; j++ )
                        {
                            for( a_inv[0][j]=0., k=j-1; k<n; k++ ) a_inv[0][j] += Q[0][k]*J_inv[k]*Q[j][k];
                            a_inv[j][0] = a_inv[0][j];
                        }
                        for( i=1; i<n; i++ )
                        {
                            for( a_inv[i][i]=0., k=i-1; k<n; k++ ) a_inv[i][i] += Q[i][k]*J_inv[k]*Q[i][k];
                            for( j=i+1; j<n; j++ )
                            {
                                for( a_inv[i][j]=0., k=j-1; k<n; k++ ) a_inv[i][j] += Q[i][k]*J_inv[k]*Q[j][k];
                                a_inv[j][i] = a_inv[i][j];
                            }
                        }
                        break;

                }//schema
                break;

            case 1: //матрица простой структуры
//                System.out.println(" simple structure matrix:");
//                System.out.println(" schema = "+schema);
                switch( schema )
                {
                    case 1:
                        //TJ
                        //первая строка
                        a[0][0] = J[0];
                        a[0][1] =-J[1];
//			for(j=2; j<n; j++ ) a[0][j] = 0.;
                        //до последней
                        for(i=1; i<n-1; i++)
                        {
//				for(j=0; j<i-1; j++ ) a[i][j] = 0.;
                            a[i][i-1] =-J[i-1];
                            a[i][i] = J[i]+J[i];
                            a[i][i+1] =-J[i+1];
//				for(j=i+2; j<n; j++ ) a[i][j] = 0.;
                        }
                        //последняя (n-1)
//			for(j=0; j<n-2; j++ ) a[n-1][j] = 0.;
                        a[n-1][n-2] =-J[n-2];
                        a[n-1][n-1] = J[n-1]+J[n-1];

                        //(TJ)T^{-1}
                        //первая строка
                        aa[1] = a[0][0];  aa[2] = a[0][1];
                        a[0][0] = aa[1]*(double)n + aa[2]*(double)(n-1);
                        double s = aa[1] + aa[2];
                        for( j=1; j<n; j++ ) a[0][j] = s*(double)(n-j);
                        //до последней
                        for(i=1; i<n-1; i++)
                        {
                            aa[0] = a[i][i-1];  aa[1] = a[i][i];  aa[2] = a[i][i+1];
                            for(j=0; j<i; j++ ) a[i][j] = aa[0]*(double)(n-i+1) + aa[1]*(double)(n-i) + aa[2]*(double)(n-i-1);
                            s = aa[0] + aa[1];
                            a[i][i] = s*(double)(n-i) + aa[2]*(double)(n-i-1);
                            s += aa[2];
                            for(j=i+1; j<n; j++ ) a[i][j] = s*(double)(n-j);
                        }
                        //последняя (n-1)
                        aa[0] = a[n-1][n-2];  aa[1] = a[n-1][n-1];
                        s = aa[0] + aa[0] + aa[1];
                        for(j=0; j<n-1; j++ ) a[n-1][j] = s;
                        a[n-1][n-1] = aa[0] + aa[1];
//_______

                        //TJ^{-1}
                        //первая строка
                        a_inv[0][0] = J_inv[0];
                        a_inv[0][1] =-J_inv[1];
                        //до последней
                        for(i=1; i<n-1; i++)
                        {
                            a_inv[i][i-1] =-J_inv[i-1];
                            a_inv[i][i] = J_inv[i]+J_inv[i];
                            a_inv[i][i+1] =-J_inv[i+1];
                        }
                        //последняя (n-1)
                        a_inv[n-1][n-2] =-J_inv[n-2];
                        a_inv[n-1][n-1] = J_inv[n-1]+J_inv[n-1];

                        //(TJ^{-1})T^{-1}
                        //первая строка
                        aa[1] = a_inv[0][0];  aa[2] = a_inv[0][1];
                        a_inv[0][0] = aa[1]*(double)n + aa[2]*(double)(n-1);
                        s = aa[1] + aa[2];
                        for( j=1; j<n; j++ ) a_inv[0][j] = s*(double)(n-j);
                        //до последней
                        for(i=1; i<n-1; i++)
                        {
                            aa[0] = a_inv[i][i-1];  aa[1] = a_inv[i][i];  aa[2] = a_inv[i][i+1];
                            for(j=0; j<i; j++ ) a_inv[i][j] = aa[0]*(double)(n-i+1) + aa[1]*(double)(n-i) + aa[2]*(double)(n-i-1);
                            s = aa[0] + aa[1];
                            a_inv[i][i] = s*(double)(n-i) + aa[2]*(double)(n-i-1);
                            s += aa[2];
                            for(j=i+1; j<n; j++ ) a_inv[i][j] = s*(double)(n-j);
                        }
                        //последняя (n-1)
                        aa[0] = a_inv[n-1][n-2];  aa[1] = a_inv[n-1][n-1];
                        s = aa[0] + aa[0] + aa[1];
                        for(j=0; j<n-1; j++ ) a_inv[n-1][j] = s;
                        a_inv[n-1][n-1] = aa[0] + aa[1];
                        break;

                }//schema
                break;

            case 2: //одна жорданова клетка 2x2 при минимальном с.з.
//                System.out.println(" J_2 type matrix: must be n > 2");
//                System.out.println(" schema = "+schema);

                switch( schema )
                {
                    case 1:
                        //TJ
                        //первая строка
                        a[0][0] = J[0];
                        a[0][1] = 1.-J[0];
                        //вторая строка
                        a[1][0] =-J[0];
                        a[1][1] =-1.+J[0]+J[0];
                        a[1][2] =-J[2];
                        //третья строка
                        a[2][1] =-J[0];
                        a[2][2] = J[2]+J[2];
                        if( n>3 ) a[2][3] =-J[3];
                        //до последней
                        for(i=3; i<n-1; i++)
                        {
                            a[i][i-1] =-J[i-1];
                            a[i][i] = J[i]+J[i];
                            a[i][i+1] =-J[i+1];
                        }
                        //последняя (n-1)
                        if( n>3 )
                        {
                            a[n-1][n-2] =-J[n-2];
                            a[n-1][n-1] = J[n-1]+J[n-1];
                        }

                        //(TJ)T^{-1}
                        //первая строка
                        aa[1] = a[0][0];  aa[2] = a[0][1];
                        a[0][0] = aa[1]*(double)n + aa[2]*(double)(n-1);
                        double s = aa[1] + aa[2];
                        for( j=1; j<n; j++ ) a[0][j] = s*(double)(n-j);
                        //до последней
                        for(i=1; i<n-1; i++)
                        {
                            aa[0] = a[i][i-1];  aa[1] = a[i][i];  aa[2] = a[i][i+1];
                            for(j=0; j<i; j++ ) a[i][j] = aa[0]*(double)(n-i+1) + aa[1]*(double)(n-i) + aa[2]*(double)(n-i-1);
                            s = aa[0] + aa[1];
                            a[i][i] = s*(double)(n-i) + aa[2]*(double)(n-i-1);
                            s += aa[2];
                            for(j=i+1; j<n; j++ ) a[i][j] = s*(double)(n-j);
                        }
                        //последняя (n-1)
                        aa[0] = a[n-1][n-2];  aa[1] = a[n-1][n-1];
                        s = aa[0] + aa[0] + aa[1];
                        for(j=0; j<n-1; j++ ) a[n-1][j] = s;
                        a[n-1][n-1] = aa[0] + aa[1];
//_______
                        //TJ^{-1}
                        //первая строка
                        a_inv[0][0] = J_inv[0];
                        a_inv[0][1] = -J_inv[0]*J_inv[0]-J_inv[0];
                        //вторая строка
                        a_inv[1][0] =-J_inv[0];
                        a_inv[1][1] = J_inv[0]*J_inv[0]+J_inv[0]+J_inv[0];
                        a_inv[1][2] =-J_inv[2];
                        //третья строка
                        a_inv[2][1] =-J_inv[0];
                        a_inv[2][2] = J_inv[2]+J_inv[2];
                        if( n>3 ) a_inv[2][3] =-J_inv[3];
                        //до последней
                        for(i=3; i<n-1; i++)
                        {
                            a_inv[i][i-1] =-J_inv[i-1];
                            a_inv[i][i] = J_inv[i]+J_inv[i];
                            a_inv[i][i+1] =-J_inv[i+1];
                        }
                        //последняя (n-1)
                        if( n>3 )
                        {
                            a_inv[n-1][n-2] =-J_inv[n-2];
                            a_inv[n-1][n-1] = J_inv[n-1]+J_inv[n-1];
                        }

                        //(TJ^{-1})T^{-1}
                        //первая строка
                        aa[1] = a_inv[0][0];  aa[2] = a_inv[0][1];
                        a_inv[0][0] = aa[1]*(double)n + aa[2]*(double)(n-1);
                        s = aa[1] + aa[2];
                        for( j=1; j<n; j++ ) a_inv[0][j] = s*(double)(n-j);
                        //до последней
                        for(i=1; i<n-1; i++)
                        {
                            aa[0] = a_inv[i][i-1];  aa[1] = a_inv[i][i];  aa[2] = a_inv[i][i+1];
                            for(j=0; j<i; j++ ) a_inv[i][j] = aa[0]*(double)(n-i+1) + aa[1]*(double)(n-i) + aa[2]*(double)(n-i-1);
                            s = aa[0] + aa[1];
                            a_inv[i][i] = s*(double)(n-i) + aa[2]*(double)(n-i-1);
                            s += aa[2];
                            for(j=i+1; j<n; j++ ) a_inv[i][j] = s*(double)(n-j);
                        }
                        //последняя (n-1)
                        aa[0] = a_inv[n-1][n-2];  aa[1] = a_inv[n-1][n-1];
                        s = aa[0] + aa[0] + aa[1];
                        for(j=0; j<n-1; j++ ) a_inv[n-1][j] = s;
                        a_inv[n-1][n-1] = aa[0] + aa[1];


                        break;
                }//schema

                break;

        }//variant

//______________________________________________________________________

        double norm,norm_inv;

        norm = matr_inf_norm ( a, n );
        genMatrixParam.setA_norm(norm);
//        System.out.println(" ||  A  || = " + norm );

        norm_inv = matr_inf_norm ( a_inv, n );
//        System.out.println(" ||A_inv|| = " + norm_inv );
        genMatrixParam.setA_inv_norm(norm_inv);
        double obusl = norm*norm_inv;
//        System.out.println(" obusl = " + obusl );
        genMatrixParam.setObusl(obusl);
        //невязка генерации
        double[][] r = matr_mul ( a, a_inv, n );

        for( i = 0; i < n; i++ ) r[i][i] -= 1.;

        norm = matr_inf_norm ( r, n );
//        System.out.println(" ||R_gen|| = " + norm );

        return genMatrixParam;
    }//mygen


}
