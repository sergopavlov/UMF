using System;
using System.IO;

namespace UMF
{
    class Program
    {
        static void Main(string[] args)
        {
            //u = t*(x-3)(x+3)
            //MKE mke = new MKE((x) => 0, (u, x) => u * u, (u, x) => 2 * u, (x, t) => t * t * x * x * x, (x) => 1, new BC1((t) => -2*t), new BC1((t) => 2 * t));//u0 sigma dersigma f lambda BCL BCR
            MKE mke = new MKE((x) => 0, (u, x) => u, (u, x) => 1, (x, t) => t * (x - 3) * (x + 3) * (x - 3) * (x + 3) - 2 * t, (x) => 1, new BC1((t) => 0), new BC1((t) => 0));//u0 sigma dersigma f lambda BCL BCR
            //MKE mke = new MKE((x) => 0, (x) => 1, (x, t) => 1, (x) => 0);
            mke.ReadMeshAndBC();
            mke.Solve(1e-15, 10000, 0.8);
            for (int i = 0; i < mke.Timegrid.Length; i++)
            {
                string res = "";
                for (int j = 0; j < mke.n; j++)
                {
                    res += $"{mke.q[i][j]} ";
                }
                Console.WriteLine(res);
            }
            Console.WriteLine("trak sosi");
        }
    }
    public class MKE
    {
        Func<double, double> u0;
        Func<double, double, double> sigma;
        Func<double, double, double> dersigma;
        Func<double, double, double> f;
        Func<double, double> lambda;
        public int n;
        private int elemcount;
        private double[] Xgrid;
        public double[] Timegrid;
        double[][] A;
        public double[][] q;
        double[] b;
        double[] Aq;
        private BC LeftBC;
        private BC RightBC;
        public MKE(Func<double, double> u0, Func<double, double, double> sigma, Func<double, double, double> dersigma, Func<double, double, double> f, Func<double, double> lambda, BC LeftBC, BC RightBC)
        {
            this.lambda = lambda;
            this.f = f;
            this.u0 = u0;
            this.sigma = sigma;
            this.dersigma = dersigma;
            this.LeftBC = LeftBC;
            this.RightBC = RightBC;
        }

        public void ReadMeshAndBC()
        {
            var str = File.ReadAllLines("Mesh/XGrid.txt");
            str = str[0].Split(' ');
            int n = str.Length;
            elemcount = n - 1;
            this.n = 2 * n - 1;
            A = new double[5][];
            for (int i = 0; i < 5; i++)
            {
                A[i] = new double[this.n];
            }
            b = new double[this.n];
            Xgrid = new double[n];
            for (int i = 0; i < n; i++)
            {
                Xgrid[i] = double.Parse(str[i]);
            }
            str = File.ReadAllLines("Mesh/TimeGrid.txt");
            str = str[0].Split(' ');
            n = str.Length;
            q = new double[n][];
            Aq = new double[this.n];
            for (int i = 0; i < n; i++)
            {
                q[i] = new double[this.n];
            }
            Timegrid = new double[n];
            q[0][0] = u0(Xgrid[0]);
            for (int i = 0; i < this.n / 2; i++)
            {
                q[0][2 * i + 1] = u0((Xgrid[i] + Xgrid[i + 1]) / 2);
                q[0][2 * i + 2] = u0(Xgrid[i + 1]);
            }
            for (int i = 0; i < n; i++)
            {
                Timegrid[i] = double.Parse(str[i]);
            }
        }
        void MakeSimpleIteration(int iter, double eps, int maxiter, double omega)
        {
            for (int i = 0; i < n; i++)
            {
                q[iter][i] = q[iter - 1][i];
            }
            int iternum = 0;
            bool flag = true;
            int[] nums = new int[3];
            double[] iks = new double[3];
            double[] funcs = new double[3];
            double dt = Timegrid[iter] - Timegrid[iter - 1];
            //итерация по неявности
            while (flag && iternum < maxiter)
            {
                //создаем матрицу
                for (int i = 0; i < 5; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        A[i][j] = 0;
                    }
                }
                for (int i = 0; i < n; i++)
                {
                    b[i] = 0;
                }
                for (int elemcounter = 0; elemcounter < elemcount; elemcounter++)
                {
                    nums[0] = 2 * elemcounter;
                    nums[1] = 2 * elemcounter + 1;
                    nums[2] = 2 * elemcounter + 2;
                    iks[0] = Xgrid[elemcounter];
                    iks[1] = (Xgrid[elemcounter] + Xgrid[elemcounter + 1]) / 2;
                    iks[2] = Xgrid[elemcounter + 1];
                    for (int i = 0; i < 3; i++)
                    {
                        funcs[i] = f(iks[i], Timegrid[iter]);
                    }
                    double avglambda = lambda((Xgrid[elemcounter] + Xgrid[elemcounter + 1]) / 2);
                    double hx = Xgrid[elemcounter + 1] - Xgrid[elemcounter];
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            A[j - i + 2][nums[i]] += avglambda * Matrices.GMatr[i][j] / hx;

                            for (int k = 0; k < 3; k++)
                            {
                                A[j - i + 2][nums[i]] += hx * sigma(Getsollution(iks[k], iter), iks[k]) * Matrices.MMatr[k][i][j] / dt;
                            }
                        }
                        double curb = 0;
                        for (int j = 0; j < 3; j++)
                        {
                            double cursigma = sigma(Getsollution(iks[j], iter), iks[j]);
                            b[nums[i]] += funcs[j] * hx * Matrices.Fmatr[i][j];
                            for (int k = 0; k < 3; k++)
                            {
                                curb += hx * cursigma * q[iter - 1][nums[k]] * Matrices.MMatr[j][k][i] / dt;
                            }
                        }
                        b[nums[i]] += curb;
                    }
                }
                //краевые условия 
                double value;
                switch (LeftBC.type)
                {
                    case 1:
                        value = ((BC1)LeftBC).value(Timegrid[iter]);
                        b[0] = value;
                        A[2][0] = 1;
                        b[1] -= value * A[3][0];
                        b[2] -= value * A[4][0];
                        A[0][2] = 0;
                        A[1][1] = 0;
                        A[3][0] = 0;
                        A[4][0] = 0;
                        break;
                    case 2:
                        value = ((BC2)LeftBC).theta(Timegrid[iter]);
                        b[0] += value;
                        break;
                    case 3:
                        value = ((BC3)LeftBC).beta(Timegrid[iter]);
                        b[0] += value * ((BC3)LeftBC).ub(Timegrid[iter]);
                        A[0][0] += value;
                        break;
                    default:
                        break;
                }
                switch (RightBC.type)
                {
                    case 1:
                        value = ((BC1)RightBC).value(Timegrid[iter]);
                        b[n - 1] = value;
                        A[2][n - 1] = 1;
                        b[n - 2] -= value * A[3][n - 2];
                        b[n - 3] -= value * A[4][n - 3];
                        A[0][n - 1] = 0;
                        A[1][n - 1] = 0;
                        A[3][n - 2] = 0;
                        A[4][n - 3] = 0;
                        break;
                    case 2://не сделано
                        //b[0] += BC1.Item2;
                        break;
                    case 3://не сделано
                        //b[0] += BC1.Item2 * BC1.Item3;
                        //A[0][0] += BC1.Item2;
                        break;
                    default:
                        break;
                }
                //проверим критерий остановки
                MatrixMult(iter);
                double summ1 = 0, summ2 = 0;
                for (int i = 0; i < n; i++)
                {
                    summ1 += (Aq[i] - b[i]) * (Aq[i] - b[i]);
                    summ2 += (b[i]) * (b[i]);
                }
                if (Math.Sqrt(summ1 / summ2) < eps)
                    flag = false;
                //решаем слау
                if (flag)
                {
                    LU();
                    Gaus(iter, omega);
                    iternum++;
                }
            }
            Console.WriteLine($"{iter} {iternum}");
        }
        void LU()
        {
            //1, 2 строчка
            A[3][0] /= A[2][0];
            for (int i = 2; i < n; i++)
            {
                double sumrow = 0;
                //нижний треугольник 
                A[1][i] -= A[0][i] * A[3][i - 2];
                //верхний
                A[4][i - 2] /= A[2][i - 2];
                A[3][i - 1] = (A[3][i - 1] - A[4][i - 2] * A[1][i - 1]) / A[2][i - 1];
                //диагональ
                A[2][i] -= A[0][i] * A[4][i - 2] + A[1][i] * A[3][i - 1];
            }
        }
        void Gaus(int iter, double omega)
        {
            b[0] /= A[2][0];
            b[1] = (b[1] - b[0] * A[1][1]) / A[2][1];
            for (int i = 2; i < n; i++)
            {
                b[i] = (b[i] - b[i - 2] * A[0][i] - b[i - 1] * A[1][i]) / A[2][i];
            }
            b[n - 2] -= A[3][n - 2] * b[n - 1];
            for (int i = n - 3; i >= 0; i--)
            {
                b[i] -= A[4][i] * b[i + 2] + A[3][i] * b[i + 1];
            }
            for (int i = 0; i < n; i++)
            {
                q[iter][i] = omega * b[i] + (1 - omega) * q[iter][i];
            }
        }
        void MatrixMult(int iter)
        {
            for (int i = 0; i < n; i++)
            {
                int jbeg = Math.Max(0, 2 - i);
                int jend = Math.Min(5, 5 + n - 3 - i);
                double cur = 0;
                for (int j = jbeg; j < jend; j++)
                {
                    cur += A[j][i] * q[iter][i - 2 + j];
                }
                Aq[i] = cur;
            }
        }
        double Getsollution(double x, int iter)
        {
            bool flag = true;
            int i = 0;
            for (; i < elemcount && flag; i++)
            {
                if (x >= Xgrid[i] && x <= Xgrid[i + 1])
                    flag = false;
            }
            i--;
            if (flag)
                throw new Exception();
            double t = (x - Xgrid[i]) / (Xgrid[i + 1] - Xgrid[i]);
            double[] psi = new double[]
            {
                2*(t-0.5)*(t-1),
                -4*t*(t-1),
                2*t*(t-0.5)
            };
            double res = 0;
            for (int p = 0; p < 3; p++)
            {
                res += q[iter][i * 2 + p] * psi[p];
            }
            return res;
        }
        public void Solve(double eps, int maxiter, double omega)
        {
            int m = Timegrid.Length;
            for (int i = 1; i < m; i++)
            {
                MakeSimpleIteration(i, eps, maxiter, omega);
            }
        }
    }
    public static class Matrices
    {
        public static double[][] GMatr = new double[][]
            {
                new double[]
                {
                    7.0/3,-8.0/3,1.0/3
                },
                new double[]
                {
                    -8.0/3,16.0/3,-8.0/3
                },
                new double[]
                {
                    1.0/3,-8.0/3,7.0/3
                }
            };
        public static double[][] Fmatr = new double[][]
        {
            new double[]
            {
                4.0/30,2.0/30,-1.0/30
            },
            new double[]
            {
                2.0/30,16.0/30,2.0/30
            },
            new double[]
            {
                -1.0/30,2.0/30,4.0/30
            }
        };
        public static double[][][] MMatr = new double[][][]
        {
            new double[][]
            {
                new double[]
                {
                    39.0/420,20.0/420,-3.0/420
                },
                new double[]
                {
                    20.0/420,16.0/420,-8.0/420
                },
                new double[]
                {
                    -3.0/420,-8.0/420,-3.0/420
                }
            },
            new double[][]
            {
                new double[]
                {
                    5.0/105,4.0/105,-2.0/105
                },
                new double[]
                {
                    4.0/105,48.0/105,4.0/105
                },
                new double[]
                {
                    -2.0/105,4.0/105,5.0/105
                }
            },
            new double[][]
            {
                new double[]
                {
                    -3.0/420,-8.0/420,-3.0/420
                },
                new double[]
                {
                    -8.0/420,16.0/420,20.0/420
                },
                new double[]
                {
                    -3.0/420,20.0/420,39.0/420
                }
            }
        };
    }
    public abstract class BC
    {
        public int type { get; init; }
        public BC(int type)
        {
            this.type = type;
        }
    }
    public class BC1 : BC
    {
        public Func<double, double> value;
        public BC1(Func<double, double> value) : base(1)
        {
            this.value = value;
        }
    }
    public class BC2 : BC
    {
        public Func<double, double> theta;
        public BC2(Func<double, double> theta) : base(2)
        {
            this.theta = theta;
        }
    }
    public class BC3 : BC
    {
        public Func<double, double> beta;
        public Func<double, double> ub;
        public BC3(Func<double, double> beta, Func<double, double> ub) : base(3)
        {
            this.beta = beta;
            this.ub = ub;
        }
    }
}
