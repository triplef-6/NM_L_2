public class Table {
    public static void printTableBeta(int n, int maxSBeta, double[] x_t, double epsilon) {
        System.out.println("""
                        +-----+-----+----------+-----------+---------+----------+-------+----------+-----+
                        !  α  !  β  !  ||A||∞  !  ||A_||∞  !  ν∞(A)  !  ||z||∞  !   ζ   !  ||R||∞  !  I  !
                        +-----+-----+----------+-----------+---------+----------+-------+----------+-----+""");
        double alpha = 1;
        double beta = 10;
        double maxBeta = Math.pow(10, maxSBeta);
        int k = 2;
        do {
            GeneratorSystem generatorSystem = new GeneratorSystem(n, x_t, alpha, beta, epsilon);

            double z_norm = z_norm(generatorSystem.getX(), generatorSystem.getX_t(), n);
            double r_norm = r_norm(generatorSystem.getA(), generatorSystem.getX(), generatorSystem.getB(), n);


            System.out.printf("|%5.0e|%5.0e|%10.4e|%11.5e|%9.3e|%10.4e|%7.1e|%9.4e|%d|",
                    alpha,
                    beta,
                    generatorSystem.getGenMatrixParam().getA_norm(),
                    generatorSystem.getGenMatrixParam().getA_inv_norm(),
                    generatorSystem.getGenMatrixParam().getObusl(),
                    z_norm,
                    z_norm / generatorSystem.getGenMatrixParam().getA_inv_norm(),
                    r_norm,
                    generatorSystem.getIter_x()
            );
            System.out.println("\n+-----+-----+----------+-----------+---------+----------+-------+----------+-----+");

            beta = Math.pow(10, k);
            k++;
        } while (beta < maxBeta); // 10 ** 12
    }

    private static double z_norm(double[] x, double[] x_t, int n) {
        double[] Z = new double[n];

        for (int i = 0; i < n; i++) {
                Z[i] = x[i] - x_t[i];
        }

        return norm(Z, n);
    }

    private static double r_norm(double[][] A, double[] x, double[] b, int n) {
        double[] Ax = multiplication_nn_n(A, x, n);
        double[] R = new double[n];

        for (int i = 0; i < n; i++) {
            R[i] = Ax[i] - b[i];
        }

        return norm(R, n);
    }

    private static double norm(double[] a, int n) {
        double norm = 0.;

        for (int j = 0; j < n; j++) {
            if (norm < Math.abs(a[j])) {
                norm = Math.abs(a[j]);
            }
        }

        return norm;
    }

    private static double[] multiplication_nn_n(double[][] A, double[] x, int n) {
        double[] res = new double[n];

        for (int i = 0; i < n; i++) {
            res[i] = 0.;
            for (int j = 0; j < n; j++) {
                res[i] += A[i][j] * x[j];
            }
        }

        return res;
    }


}
