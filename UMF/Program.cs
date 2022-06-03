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
            Func<double, double, double> u = (double x, double t) => x * t;
            MKE mke = new MKE((x) => x, (u, x) => u * u, (u, x) => 2 * u, (x, t) => x * t * t * x * x, (x) => 1, new BC1((t) => t * 1), new BC1((t) => 20 * t));//u0 sigma dersigma f lambda BCL BCR
            //MKE mke = new MKE((x) => 0, (x) => 1, (x, t) => 1, (x) => 0);
            mke.ReadMesh();
            mke.SolveNewton(1e-15, 1000, 0.6);
            //mke.SolveSimpleIteration(1e-15, 10000, 0.5);
            Console.WriteLine("hello world");
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
        double[][] AL;
        double[][] A;
        public double[][] q;
        double[] bL;
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

        public void ReadMesh()
        {
            var str = File.ReadAllLines("Mesh/XGrid.txt");
            str = str[0].Split(' ');
            int n = str.Length;
            elemcount = n - 1;
            this.n = 2 * n - 1;
            AL = new double[5][];
            A = new double[5][];
            for (int i = 0; i < 5; i++)
            {
                AL[i] = new double[this.n];
                A[i] = new double[this.n];
            }
            bL = new double[this.n];
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
            double summ1 = 0, summ2 = 0;
            //итерация по неявности
            double relresidual = 1;
            while (relresidual > eps && iternum < maxiter)
            {
                //создаем матрицу
                for (int i = 0; i < 5; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        AL[i][j] = 0;
                        A[i][j] = 0;
                    }
                }
                for (int i = 0; i < n; i++)
                {
                    bL[i] = 0;
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
                            AL[j - i + 2][nums[i]] += avglambda * Matrices.GMatr[i][j] / hx;
                            A[j - i + 2][nums[i]] += avglambda * Matrices.GMatr[i][j] / hx;

                            for (int k = 0; k < 3; k++)
                            {
                                AL[j - i + 2][nums[i]] += hx * sigma(Getsollution(iks[k], iter), iks[k]) * Matrices.MMatr[k][i][j] / dt;
                                A[j - i + 2][nums[i]] += hx * sigma(Getsollution(iks[k], iter), iks[k]) * Matrices.MMatr[k][i][j] / dt;
                            }
                        }
                        double curb = 0;
                        for (int j = 0; j < 3; j++)
                        {
                            double cursigma = sigma(Getsollution(iks[j], iter), iks[j]);
                            bL[nums[i]] += funcs[j] * hx * Matrices.Fmatr[i][j];
                            b[nums[i]] += funcs[j] * hx * Matrices.Fmatr[i][j];
                            for (int k = 0; k < 3; k++)
                            {
                                curb += hx * cursigma * q[iter - 1][nums[k]] * Matrices.MMatr[j][k][i] / dt;
                            }
                        }
                        bL[nums[i]] += curb;
                        b[nums[i]] += curb;
                    }
                }
                //краевые условия 
                double value;
                switch (LeftBC.type)
                {
                    case 1:
                        value = ((BC1)LeftBC).value(Timegrid[iter]);
                        bL[0] = value;
                        AL[2][0] = 1;
                        bL[1] -= value * AL[3][0];
                        bL[2] -= value * AL[4][0];
                        AL[0][2] = 0;
                        AL[1][1] = 0;
                        AL[3][0] = 0;
                        AL[4][0] = 0;
                        b[0] = value;
                        A[2][0] = 1;
                        b[1] -= value * AL[3][0];
                        b[2] -= value * AL[4][0];
                        A[0][2] = 0;
                        A[1][1] = 0;
                        A[3][0] = 0;
                        A[4][0] = 0;
                        break;
                    case 2:
                        value = ((BC2)LeftBC).theta(Timegrid[iter]);
                        bL[0] += value;
                        break;
                    case 3:
                        value = ((BC3)LeftBC).beta(Timegrid[iter]);
                        bL[0] += value * ((BC3)LeftBC).ub(Timegrid[iter]);
                        AL[0][0] += value;
                        break;
                    default:
                        break;
                }
                switch (RightBC.type)
                {
                    case 1:
                        value = ((BC1)RightBC).value(Timegrid[iter]);
                        bL[n - 1] = value;
                        AL[2][n - 1] = 1;
                        bL[n - 2] -= value * AL[3][n - 2];
                        bL[n - 3] -= value * AL[4][n - 3];
                        AL[0][n - 1] = 0;
                        AL[1][n - 1] = 0;
                        AL[3][n - 2] = 0;
                        AL[4][n - 3] = 0;
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
                //решаем слау
                MatrixMult(iter);
                summ1 = 0;
                summ2 = 0;
                for (int i = 0; i < n; i++)
                {
                    summ1 += (Aq[i] - bL[i]) * (Aq[i] - bL[i]);
                    summ2 += (bL[i]) * (bL[i]);
                }
                relresidual = Math.Sqrt(summ1 / summ2);
                LU();
                Gaus(iter, omega);
                iternum++;

            }
            Console.WriteLine($"{iter} {iternum} {Math.Sqrt(summ1 / summ2)}");
        }
        void MakeNewtonIteration(int iter, double eps, int maxiter, double omega)
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
            double relresidual = 1;
            while (relresidual > eps && iternum < maxiter)
            {
                //создаем матрицу
                for (int i = 0; i < 5; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        AL[i][j] = 0;
                        A[i][j] = 0;
                    }
                }
                for (int i = 0; i < n; i++)
                {
                    bL[i] = 0;
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
                            AL[j - i + 2][nums[i]] += avglambda * Matrices.GMatr[i][j] / hx;
                            A[j - i + 2][nums[i]] += avglambda * Matrices.GMatr[i][j] / hx;
                            for (int k = 0; k < 3; k++)
                            {
                                AL[j - i + 2][nums[i]] += hx * sigma(Getsollution(iks[k], iter), iks[k]) * Matrices.MMatr[k][i][j] / dt;//вклад от A(q0)
                                A[j - i + 2][nums[i]] += hx * sigma(Getsollution(iks[k], iter), iks[k]) * Matrices.MMatr[k][i][j] / dt;//вклад от A(q0)

                            }
                            for (int r = 0; r < 3; r++)
                            {
                                for (int p = 0; p < 3; p++)
                                {

                                    AL[j - i + 2][nums[i]] += hx / dt * dersigma(Getsollution(iks[p], iter), iks[p]) * Matrices.MMatr2[j][i][r][p] * q[iter][nums[r]];
                                }
                                //AL[j - i + 2][nums[i]] += hx / dt * dersigma(q[iter][nums[j]], iks[j]) * Matrices.MMatr[j][i][r] * q[iter][nums[r]];
                            }
                            for (int k = 0; k < 3; k++)//вклад от линеаризации правой части
                            {
                                for (int p = 0; p < 3; p++)
                                {
                                    AL[j - i + 2][nums[i]] -= hx / dt * dersigma(Getsollution(iks[p], iter), iks[p]) * Matrices.MMatr2[p][k][i][j] * q[iter - 1][nums[k]];

                                }
                                //AL[j - i + 2][nums[i]] -= hx / dt * dersigma(q[iter][nums[j]], iks[j]) * Matrices.MMatr[k][i][j] * q[iter - 1][nums[k]];
                            }
                        }
                        double curb = 0;
                        for (int j = 0; j < 3; j++)
                        {
                            double cursigma = sigma(Getsollution(iks[j], iter), iks[j]);
                            bL[nums[i]] += funcs[j] * hx * Matrices.Fmatr[i][j];
                            b[nums[i]] += funcs[j] * hx * Matrices.Fmatr[i][j];

                            for (int k = 0; k < 3; k++)
                            {
                                curb += hx * cursigma * q[iter - 1][nums[k]] * Matrices.MMatr[j][k][i] / dt;
                            }
                        }
                        bL[nums[i]] += curb;
                        b[nums[i]] += curb;
                        for (int j = 0; j < 3; j++)
                        {
                            for (int r = 0; r < 3; r++)
                            {
                                for (int p = 0; p < 3; p++)
                                {
                                    bL[nums[i]] += hx / dt * q[iter][nums[r]] * q[iter][nums[j]] * dersigma(q[iter][nums[p]], iks[p]) * Matrices.MMatr2[p][i][j][r];

                                }
                                //bL[nums[i]] += hx / dt * q[iter][nums[r]] * q[iter][nums[j]] * dersigma(q[iter][nums[r]], iks[r]) * Matrices.MMatr[i][j][r];
                            }
                        }
                        for (int r = 0; r < 3; r++)
                        {
                            for (int k = 0; k < 3; k++)//цикл умноженич на матрицу
                            {
                                for (int p = 0; p < 3; p++)
                                {
                                    bL[nums[i]] -= hx / dt * q[iter - 1][nums[k]] * q[iter][nums[r]] * dersigma(q[iter][nums[p]], iks[p]) * Matrices.MMatr2[p][i][k][r];

                                }
                                //bL[nums[i]] -= hx / dt * q[iter - 1][nums[k]] * q[iter][nums[r]] * dersigma(q[iter][nums[r]], iks[r]) * Matrices.MMatr[i][k][r];

                            }
                        }

                    }
                }
                //краевые условия 
                double value;
                switch (LeftBC.type)
                {
                    case 1:
                        value = ((BC1)LeftBC).value(Timegrid[iter]);
                        bL[0] = value;
                        AL[2][0] = 1;
                        bL[1] -= value * AL[3][0];
                        bL[2] -= value * AL[4][0];
                        AL[0][2] = 0;
                        AL[1][1] = 0;
                        AL[3][0] = 0;
                        AL[4][0] = 0;
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
                        bL[0] += value;
                        break;
                    case 3:
                        value = ((BC3)LeftBC).beta(Timegrid[iter]);
                        bL[0] += value * ((BC3)LeftBC).ub(Timegrid[iter]);
                        AL[0][0] += value;
                        break;
                    default:
                        break;
                }
                switch (RightBC.type)
                {
                    case 1:
                        value = ((BC1)RightBC).value(Timegrid[iter]);
                        bL[n - 1] = value;
                        AL[2][n - 1] = 1;
                        bL[n - 2] -= value * AL[3][n - 2];
                        bL[n - 3] -= value * AL[4][n - 3];
                        AL[0][n - 1] = 0;
                        AL[1][n - 1] = 0;
                        AL[3][n - 2] = 0;
                        AL[4][n - 3] = 0;
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
                relresidual = Math.Sqrt(summ1 / summ2);

                //решаем слау
                LU();
                Gaus(iter, omega);
                iternum++;
            }
            Console.WriteLine($"{iter} {iternum} {relresidual}");
        }
        void LU()
        {
            //1, 2 строчка
            AL[3][0] /= AL[2][0];
            for (int i = 2; i < n; i++)
            {
                //нижний треугольник 
                AL[1][i] -= AL[0][i] * AL[3][i - 2];
                //верхний
                AL[4][i - 2] /= AL[2][i - 2];
                AL[3][i - 1] = (AL[3][i - 1] - AL[4][i - 2] * AL[1][i - 1]) / AL[2][i - 1];
                //диагональ
                AL[2][i] -= AL[0][i] * AL[4][i - 2] + AL[1][i] * AL[3][i - 1];
            }
        }
        void Gaus(int iter, double omega)
        {
            bL[0] /= AL[2][0];
            bL[1] = (bL[1] - bL[0] * AL[1][1]) / AL[2][1];
            for (int i = 2; i < n; i++)
            {
                bL[i] = (bL[i] - bL[i - 2] * AL[0][i] - bL[i - 1] * AL[1][i]) / AL[2][i];
            }
            bL[n - 2] -= AL[3][n - 2] * bL[n - 1];
            for (int i = n - 3; i >= 0; i--)
            {
                bL[i] -= AL[4][i] * bL[i + 2] + AL[3][i] * bL[i + 1];
            }
            for (int i = 0; i < n; i++)
            {
                q[iter][i] = omega * bL[i] + (1 - omega) * q[iter][i];
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
        public void SolveSimpleIteration(double eps, int maxiter, double omega)
        {
            int m = Timegrid.Length;
            for (int i = 1; i < m; i++)
            {
                MakeSimpleIteration(i, eps, maxiter, omega);
            }
        }
        public void SolveNewton(double eps, int maxiter, double omega)
        {
            int m = Timegrid.Length;
            for (int i = 1; i < m; i++)
            {
                MakeNewtonIteration(i, eps, maxiter, omega);
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
        public static double[][][][] MMatr2 = new double[][][][]
        {
            new double[][][]
{
new double[][]
{
new double[]
{
0.07301587301587302, 0.025396825396825404, -0.005555555555555557,
},
new double[]
{
0.025396825396825404, 0.0253968253968254, -0.003174603174603175,
},
new double[]
{
-0.005555555555555557, -0.0031746031746031746, 0.0015873015873015877,
},
},
new double[][]
{
new double[]
{
0.025396825396825404, 0.0253968253968254, -0.0031746031746031746,
},
new double[]
{
0.0253968253968254, 0.025396825396825383, -0.0126984126984127,
},
new double[]
{
-0.0031746031746031746, -0.0126984126984127, -0.0031746031746031755,
},
},
new double[][]
{
new double[]
{
-0.005555555555555557, -0.0031746031746031746, 0.001587301587301588,
},
new double[]
{
-0.0031746031746031746, -0.0126984126984127, -0.0031746031746031755,
},
new double[]
{
0.0015873015873015877, -0.0031746031746031755, -0.005555555555555557,
},
},
},
new double[][][]
{
new double[][]
{
new double[]
{
0.025396825396825404, 0.0253968253968254, -0.0031746031746031746,
},
new double[]
{
0.0253968253968254, 0.025396825396825383, -0.0126984126984127,
},
new double[]
{
-0.0031746031746031746, -0.0126984126984127, -0.0031746031746031755,
},
},
new double[][]
{
new double[]
{
0.025396825396825397, 0.025396825396825376, -0.0126984126984127,
},
new double[]
{
0.025396825396825383, 0.40634920634920635, 0.025396825396825414,
},
new double[]
{
-0.0126984126984127, 0.025396825396825418, 0.025396825396825404,
},
},
new double[][]
{
new double[]
{
-0.0031746031746031746, -0.0126984126984127, -0.0031746031746031755,
},
new double[]
{
-0.0126984126984127, 0.025396825396825425, 0.025396825396825407,
},
new double[]
{
-0.0031746031746031755, 0.025396825396825407, 0.025396825396825404,
},
},
},
new double[][][]
{
new double[][]
{
new double[]
{
-0.005555555555555557, -0.0031746031746031746, 0.001587301587301588,
},
new double[]
{
-0.0031746031746031746, -0.0126984126984127, -0.0031746031746031755,
},
new double[]
{
0.0015873015873015877, -0.0031746031746031755, -0.005555555555555557,
},
},
new double[][]
{
new double[]
{
-0.0031746031746031746, -0.0126984126984127, -0.0031746031746031755,
},
new double[]
{
-0.0126984126984127, 0.025396825396825425, 0.025396825396825407,
},
new double[]
{
-0.0031746031746031755, 0.025396825396825407, 0.025396825396825404,
},
},
new double[][]
{
new double[]
{
0.0015873015873015877, -0.0031746031746031755, -0.005555555555555557,
},
new double[]
{
-0.0031746031746031755, 0.025396825396825407, 0.025396825396825404,
},
new double[]
{
-0.0055555555555555575, 0.025396825396825404, 0.07301587301587298,
},
},
},
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
