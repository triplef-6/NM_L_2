import java.util.ArrayList;
import java.util.List;

public class GeneratorSystem {
    private int n;
    private double[][] A;
    private double[] b;
    private double[] x_t;
    private double[] x;
    private GenMatrixParam genMatrixParam;
    private int iter_x;

    public GeneratorSystem(int n, double[] x_t, double alpha, double beta, double epsilon) {
        Gen g = new Gen();
        this.n = n;
        this.x_t = x_t;

        A = new double[n][n];
        double[][] _A_ = new double[n][n];
        genMatrixParam = g.mygen(A, _A_, n, alpha, beta, 1, 2, 0, 1, false);

        b = multiplication_nn_n(A, x_t);
        x = system_solution(epsilon);
    }

    public double[] system_solution(double epsilon) { // сопряжённых градиентов
        double[] x0 = new double[n];
        for (double x_from_x0 : x0) {
            x_from_x0 = 0.;
        }

        double[] r0 = subtraction_n(b, multiplication_nn_n(A, x0));

        List<double[]> Xs = new ArrayList<>();
        Xs.add(x0);
        List<double[]> Rs = new ArrayList<>();
        Rs.add(r0);
        List<double[]> Zs = new ArrayList<>();
        Zs.add(r0);

        int k = 0;
        double norma_b = norma_n(b);
        do {
            k++;

            double alfaK = sca_multiplication(Rs.get(k - 1), Rs.get(k - 1)) /
                    sca_multiplication(multiplication_nn_n(A, Zs.get(k - 1)), Zs.get(k - 1));
            Xs.add(addition_n(Xs.get(k - 1), multiplication_n_a(Zs.get(k - 1), alfaK)));
            Rs.add(subtraction_n(Rs.get(k - 1), multiplication_n_a(multiplication_nn_n(A, Zs.get(k - 1)), alfaK)));
            double beta = sca_multiplication(Rs.get(k), Rs.get(k)) /
                    sca_multiplication(Rs.get(k - 1), Rs.get(k - 1));
            Zs.add(addition_n(Rs.get(k), multiplication_n_a(Zs.get(k - 1), beta)));


        } while (norma_n(Rs.get(k)) / norma_b >= epsilon);
        iter_x = k;
        return Xs.get(k);
    }

    private double[] multiplication_nn_n(double[][] A, double[] x) {
        double[] res = new double[n];

        for (int i = 0; i < n; i++) {
            res[i] = 0.;
            for (int j = 0; j < n; j++) {
                res[i] += A[i][j] * x[j];
            }
        }

        return res;
    }

    private double[] multiplication_n_a(double[] a, double b) {
        double[] res = new double[n];

        for (int i = 0; i < n; i++) {
            res[i] = a[i] * b;
        }

        return res;
    }

    private double[] addition_n(double[] a, double[] b) {
        double[] res = new double[n];

        for (int i = 0; i < n; i++) {
            res[i] = a[i] + b[i];
        }

        return res;
    }

    private double[] subtraction_n(double[] a, double[] b) {
        double[] res = new double[n];

        for (int i = 0; i < n; i++) {
            res[i] = a[i] - b[i];
        }

        return res;
    }

    private double norma_n(double[] a) {
        double norm = 0.;

        for (int i = 0; i < n; i++) {
            double s = Math.abs(a[i]);
            if (s > norm) {
                norm = s;
            }
        }

        return norm;
    }

    private double sca_multiplication(double[] a, double[] b) {
        double res = 0;

        for (int i = 0; i < n; i++) {
            res += a[i] * b[i];
        }

        return res;
    }
    public double[][] getA() {
        return A;
    }

    public double[] getB() {
        return b;
    }

    public double[] getX() {
        return x;
    }

    public double[] getX_t() {
        return x_t;
    }

    public GenMatrixParam getGenMatrixParam() {
        return genMatrixParam;
    }

    public int getN() {
        return n;
    }

    public int getIter_x() {
        return iter_x;
    }
}









