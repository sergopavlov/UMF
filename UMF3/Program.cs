using System;
using System.Collections.Generic;
using System.IO;

namespace UMF3
{
    public enum edges
    {
        Top,
        Bottom,
        Left,
        Right,
        Front,
        Back
    }
    internal class Program
    {
        static void Main(string[] args)
        {

            #region mke
            //Func<double, double, double, double> ressin = (x, y, z) => y;
            //Func<double, double, double, double> rescos = (x, y, z) => x;
            //List<BC> BCs = new();
            //for (int i = 0; i < 6; i++)
            //{
            //    BCs.Add(new BC1((edges)i, rescos, 0));
            //    BCs.Add(new BC1((edges)i, ressin, 1));
            //}
            //double omega = 10;
            //double lambda = 100000;
            //double sigma = 1;
            //double hi = 1e-11;
            //Mke mke = new Mke((x, y, z) => lambda, (x, y, z) => sigma, (x, y, z) => hi, (x, y, z) => -hi * omega * omega * x + omega * sigma * y, (x, y, z) => -hi * omega * omega * y - omega * sigma * x, omega, BCs);
            //mke.ReadMesh();
            ////mke.SolveLU();
            ////mke.SolveLOS(1e-15,10000);
            //// mke.SolveLOSPrecond(1e-15, 10000);
            //mke.SolveGMRES(1e-15, 1000000, 5);
            ////MatrixSparce test = new();
            ////test.al = new List<double>() { 3, 3, 4, 13, 3, 5, 16, 5, 7, 3 };
            ////test.au = new List<double>() { 2, 12, 5, 14, 8, 9, 5, 19, 34, 52 };
            ////test.di = new List<double>() { 1, 4, 6, 4, 3 };
            ////test.b = new List<double>() { 122, 153, 235, 310, 74 };
            ////test.ia = new List<int>() { 0, 0, 1, 3, 6, 10 };
            ////test.ja = new List<int>() { 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 };
            ////test.n = 5;
            ////List<double> x0 = new List<double>();
            ////for (int i = 0; i < 5; i++)
            ////{
            ////    x0.Add(0);
            ////}
            //double pogr = 0;
            //double norm = 0;
            //double z = 0;
            //for (int i = 0; i < mke.Zgrid.Count; i++)
            //{
            //    double y = 0;
            //    for (int j = 0; j < mke.Ygrid.Count; j++)
            //    {
            //        double x = 0;
            //        for (int k = 0; k < mke.Xgrid.Count; k++)
            //        {
            //            int index = 2 * (k + mke.Xgrid.Count * j + mke.Xgrid.Count * mke.Ygrid.Count * i);
            //            pogr += (mke.q[index] - rescos(mke.Xgrid[k], mke.Ygrid[j], mke.Zgrid[i])) * (mke.q[index] - rescos(mke.Xgrid[k], mke.Ygrid[j], mke.Zgrid[i]));
            //            pogr += (mke.q[index + 1] - ressin(mke.Xgrid[k], mke.Ygrid[j], mke.Zgrid[i])) * (mke.q[index + 1] - ressin(mke.Xgrid[k], mke.Ygrid[j], mke.Zgrid[i]));
            //            norm += rescos(mke.Xgrid[k], mke.Ygrid[j], mke.Zgrid[i]) * rescos(mke.Xgrid[k], mke.Ygrid[j], mke.Zgrid[i]);
            //            norm += ressin(mke.Xgrid[k], mke.Ygrid[j], mke.Zgrid[i]) * ressin(mke.Xgrid[k], mke.Ygrid[j], mke.Zgrid[i]);
            //            x += 1;
            //        }
            //        y += 1;
            //    }
            //    z += 1;
            //}
            //Console.WriteLine($"{Math.Sqrt(pogr / norm)} {mke.q.Count} {lambda} {sigma} {hi} {omega}");
            #endregion
        }
    }
    public class Mke
    {
        int n;
        int m;
        public List<double> Xgrid = new();
        public List<double> Ygrid = new();
        public List<double> Zgrid = new();
        MatrixSparce mat = new();
        MatrixProfile matprof;
        public List<double> q = new();
        private Func<double, double, double, double> lambda;
        private Func<double, double, double, double> sigma;
        private Func<double, double, double, double> hi;
        private Func<double, double, double, double> fc;
        private Func<double, double, double, double> fs;
        public List<Elem> elems = new();
        private List<BC> BCs;
        private double omega;
        private int GetGlobalNum(Elem el, int localnum)
        {
            return 2 * el.xlownum + 2 * n * el.ylownum + 2 * m * el.zlownum + localnum % 2
                + (localnum / 2) % 2 * 2 + (localnum / 4) % 2 * n * 2 + localnum / 8 * m * 2;
        }
        public Mke(Func<double, double, double, double> lambda, Func<double, double, double, double> sigma, Func<double, double, double, double> hi, Func<double, double, double, double> fc, Func<double, double, double, double> fs, double omega, List<BC> BCs)
        {
            this.lambda = lambda;
            this.sigma = sigma;
            this.hi = hi;
            this.fc = fc;
            this.fs = fs;
            this.omega = omega;
            this.BCs = BCs;
        }
        public void ReadMesh()
        {
            foreach (var item in File.ReadAllLines("TXT/Xgrid.txt"))
            {
                Xgrid.Add(double.Parse(item));
            }
            foreach (var item in File.ReadAllLines("TXT/Ygrid.txt"))
            {
                Ygrid.Add(double.Parse(item));
            }
            foreach (var item in File.ReadAllLines("TXT/Zgrid.txt"))
            {
                Zgrid.Add(double.Parse(item));
            }
            n = Xgrid.Count;
            m = Xgrid.Count * Ygrid.Count;
            for (int i = 0; i < n - 1; i++)
            {
                for (int j = 0; j < Ygrid.Count - 1; j++)
                {
                    for (int p = 0; p < Zgrid.Count - 1; p++)
                    {
                        elems.Add(new Elem(i, j, p));
                    }
                }
            }
        }
        public void GenerateProfile()
        {
            mat.n = 2 * Xgrid.Count * Ygrid.Count * Zgrid.Count;
            mat.di = new List<double>(mat.n);
            mat.ia = new List<int>(mat.n + 1);
            mat.ja = new List<int>();
            mat.al = new List<double>();
            mat.au = new List<double>();
            mat.b = new List<double>(mat.n);
            List<SortedSet<int>> list = new List<SortedSet<int>>(mat.n);
            for (int i = 0; i < mat.n; i++)
            {
                list.Add(new SortedSet<int>());
            }
            foreach (var item in elems)
            {
                for (int i = 0; i < 16; i++)
                {
                    for (int j = 0; j < 16; j++)
                    {
                        list[GetGlobalNum(item, i)].Add(GetGlobalNum(item, j));
                    }
                }
            }
            mat.ia.Add(0);
            mat.ia.Add(0);
            mat.di.Add(0);
            mat.b.Add(0);
            for (int i = 1; i < mat.n; i++)
            {
                mat.di.Add(0);
                mat.b.Add(0);
                int count = 0;
                foreach (var item in list[i])
                {
                    if (item < i)
                    {
                        mat.ja.Add(item);
                        mat.al.Add(0);
                        mat.au.Add(0);
                        count++;
                    }
                }
                mat.ia.Add(mat.ia[i] + count);
            }
        }
        public void Addlocal(Elem el)
        {
            double hx = Xgrid[el.xlownum + 1] - Xgrid[el.xlownum];
            double hy = Ygrid[el.ylownum + 1] - Ygrid[el.ylownum];
            double hz = Zgrid[el.zlownum + 1] - Zgrid[el.zlownum];
            double Xmean = (Xgrid[el.xlownum + 1] + Xgrid[el.xlownum]) / 2;
            double Ymean = (Ygrid[el.ylownum + 1] + Ygrid[el.ylownum]) / 2;
            double Zmean = (Zgrid[el.zlownum + 1] + Zgrid[el.zlownum]) / 2;
            double[] fcos = new double[8];
            double[] fsin = new double[8];
            for (int i = 0; i < 8; i++)
            {
                fcos[i] = fc(Xgrid[el.xlownum + i % 2], Ygrid[el.ylownum + (i / 2) % 2], Zgrid[el.zlownum + i / 4]);
                fsin[i] = fs(Xgrid[el.xlownum + i % 2], Ygrid[el.ylownum + (i / 2) % 2], Zgrid[el.zlownum + i / 4]);
            }
            for (int i = 0; i < 16; i += 2)//заполняем косинусовые строчки
            {
                double curb = 0;
                for (int j = 0; j < 8; j++)
                {
                    curb += fcos[j] * Matrices.Mij(i / 2, j, hx, hy, hz);
                }
                mat.b[GetGlobalNum(el, i)] += curb;
                mat.di[GetGlobalNum(el, i)] += lambda(Xmean, Ymean, Zmean) * Matrices.Gij(i / 2, i / 2, hx, hy, hz) - omega * omega * hi(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, i / 2, hx, hy, hz);
                for (int j = 0; j < i; j += 2)//заполняем cos cos
                {
                    int max = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, i) : GetGlobalNum(el, j);
                    int min = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, j) : GetGlobalNum(el, i);
                    int k = mat.ia[max];
                    while (mat.ja[k] < min)
                    {
                        k++;
                    }
                    mat.al[k] += lambda(Xmean, Ymean, Zmean) * Matrices.Gij(i / 2, j / 2, hx, hy, hz) - omega * omega * hi(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                    mat.au[k] += lambda(Xmean, Ymean, Zmean) * Matrices.Gij(i / 2, j / 2, hx, hy, hz) - omega * omega * hi(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                }
                for (int j = 1; j < i; j += 2)//заполняем cos sin
                {
                    int max = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, i) : GetGlobalNum(el, j);
                    int min = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, j) : GetGlobalNum(el, i);
                    int k = mat.ia[max];
                    while (mat.ja[k] < min)
                    {
                        k++;
                    }
                    mat.al[k] += omega * sigma(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                    mat.au[k] += -omega * sigma(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                }
            }
            for (int i = 1; i < 16; i += 2)//заполняем синусовые строчки
            {
                double curb = 0;
                for (int j = 0; j < 8; j++)
                {
                    curb += fsin[j] * Matrices.Mij(i / 2, j, hx, hy, hz);
                }
                mat.b[GetGlobalNum(el, i)] += curb;
                mat.di[GetGlobalNum(el, i)] += lambda(Xmean, Ymean, Zmean) * Matrices.Gij(i / 2, i / 2, hx, hy, hz) - omega * omega * hi(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, i / 2, hx, hy, hz);
                for (int j = 0; j < i; j += 2)// заполняем sin cos
                {
                    int max = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, i) : GetGlobalNum(el, j);
                    int min = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, j) : GetGlobalNum(el, i);
                    int k = mat.ia[max];
                    while (mat.ja[k] < min)
                    {
                        k++;
                    }
                    mat.al[k] += -omega * sigma(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                    mat.au[k] += omega * sigma(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                }
                for (int j = 1; j < i; j += 2)// заполняем sin sin
                {
                    int max = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, i) : GetGlobalNum(el, j);
                    int min = GetGlobalNum(el, i) > GetGlobalNum(el, j) ? GetGlobalNum(el, j) : GetGlobalNum(el, i);
                    int k = mat.ia[max];
                    while (mat.ja[k] < min)
                    {
                        k++;
                    }
                    mat.al[k] += lambda(Xmean, Ymean, Zmean) * Matrices.Gij(i / 2, j / 2, hx, hy, hz) - omega * omega * hi(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                    mat.au[k] += lambda(Xmean, Ymean, Zmean) * Matrices.Gij(i / 2, j / 2, hx, hy, hz) - omega * omega * hi(Xmean, Ymean, Zmean) * Matrices.Mij(i / 2, j / 2, hx, hy, hz);
                }
            }
        }
        public void AddBoundary3()
        {
            foreach (var item in BCs.FindAll((T) => T.bctype == 3))
            {
                var curitem = item as BC3;
                switch (curitem.edge)
                {
                    case edges.Top:
                        foreach (var el in elems.FindAll((T) => T.zlownum + 1 == Zgrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 8), GetGlobalNum(el, 10), GetGlobalNum(el, 12), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 9), GetGlobalNum(el, 11), GetGlobalNum(el, 13), GetGlobalNum(el, 15) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double z = Zgrid[Zgrid.Count - 1];
                            double hx = x[1] - x[0], hy = y[1] - y[0];

                            for (int i = 0; i < 4; i++)
                            {
                                double locb = 0;
                                for (int j = 0; j < 4; j++)
                                {
                                    locb += curitem.betta * curitem.value(x[j % 2], y[j / 2], z) * Matrices.Mij2(i, j, hx, hy);
                                }
                                mat.b[index[i]] += locb;
                                mat.di[index[i]] += curitem.betta * Matrices.Mij2(i, i, hx, hy);
                                for (int j = 0; j < i; j++)
                                {
                                    int k = mat.ia[index[i]];
                                    while (mat.ja[k] < index[j])
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[j])
                                    {
                                        mat.al[k] += curitem.betta * Matrices.Mij2(i, j, hx, hy);
                                        mat.au[k] += curitem.betta * Matrices.Mij2(i, j, hx, hy);
                                    }
                                }
                            }
                        }
                        break;
                    case edges.Bottom:
                        foreach (var el in elems.FindAll((T) => T.zlownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 2), GetGlobalNum(el, 4), GetGlobalNum(el, 6) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 3), GetGlobalNum(el, 5), GetGlobalNum(el, 7) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double z = Zgrid[0];
                            double hx = x[1] - x[0], hy = y[1] - y[0];

                            for (int i = 0; i < 4; i++)
                            {
                                double locb = 0;
                                for (int j = 0; j < 4; j++)
                                {
                                    locb += curitem.betta * curitem.value(x[j % 2], y[j / 2], z) * Matrices.Mij2(i, j, hx, hy);
                                }
                                mat.b[index[i]] += locb;
                                mat.di[index[i]] += curitem.betta * Matrices.Mij2(i, i, hx, hy);
                                for (int j = 0; j < i; j++)
                                {
                                    int k = mat.ia[index[i]];
                                    while (mat.ja[k] < index[j])
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[j])
                                    {
                                        mat.al[k] += curitem.betta * Matrices.Mij2(i, j, hx, hy);
                                        mat.au[k] += curitem.betta * Matrices.Mij2(i, j, hx, hy);
                                    }
                                }
                            }
                        }
                        break;
                    case edges.Left:
                        foreach (var el in elems.FindAll((T) => T.xlownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 4), GetGlobalNum(el, 8), GetGlobalNum(el, 12) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 5), GetGlobalNum(el, 9), GetGlobalNum(el, 13) };
                            }
                            double x = Xgrid[0];
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hy = y[1] - y[0];

                            for (int i = 0; i < 4; i++)
                            {
                                double locb = 0;
                                for (int j = 0; j < 4; j++)
                                {
                                    locb += curitem.betta * curitem.value(x, y[j % 2], z[j / 2]) * Matrices.Mij2(i, j, hz, hy);
                                }
                                mat.b[index[i]] += locb;
                                mat.di[index[i]] += curitem.betta * Matrices.Mij2(i, i, hz, hy);
                                for (int j = 0; j < i; j++)
                                {
                                    int k = mat.ia[index[i]];
                                    while (mat.ja[k] < index[j])
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[j])
                                    {
                                        mat.al[k] += curitem.betta * Matrices.Mij2(i, j, hz, hy);
                                        mat.au[k] += curitem.betta * Matrices.Mij2(i, j, hz, hy);
                                    }
                                }
                            }


                        }
                        break;
                    case edges.Right:
                        foreach (var el in elems.FindAll((T) => T.xlownum + 1 == Xgrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 2), GetGlobalNum(el, 6), GetGlobalNum(el, 10), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 3), GetGlobalNum(el, 7), GetGlobalNum(el, 11), GetGlobalNum(el, 15) };
                            }
                            double x = Xgrid[Xgrid.Count - 1];
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hy = y[1] - y[0];
                            for (int i = 0; i < 4; i++)
                            {
                                double locb = 0;
                                for (int j = 0; j < 4; j++)
                                {
                                    locb += curitem.betta * curitem.value(x, y[j % 2], z[j / 2]) * Matrices.Mij2(i, j, hz, hy);
                                }
                                mat.b[index[i]] += locb;
                                mat.di[index[i]] += curitem.betta * Matrices.Mij2(i, i, hz, hy);
                                for (int j = 0; j < i; j++)
                                {
                                    int k = mat.ia[index[i]];
                                    while (mat.ja[k] < index[j])
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[j])
                                    {
                                        mat.al[k] += curitem.betta * Matrices.Mij2(i, j, hz, hy);
                                        mat.au[k] += curitem.betta * Matrices.Mij2(i, j, hz, hy);
                                    }
                                }
                            }
                        }
                        break;
                    case edges.Front:
                        foreach (var el in elems.FindAll((T) => T.ylownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 2), GetGlobalNum(el, 8), GetGlobalNum(el, 10) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 3), GetGlobalNum(el, 9), GetGlobalNum(el, 11) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double y = Ygrid[0];
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hx = x[1] - x[0];

                            for (int i = 0; i < 4; i++)
                            {
                                double locb = 0;
                                for (int j = 0; j < 4; j++)
                                {
                                    locb += curitem.betta * curitem.value(x[j % 2], y, z[j / 2]) * Matrices.Mij2(i, j, hz, hz);
                                }
                                mat.b[index[i]] += locb;
                                mat.di[index[i]] += curitem.betta * Matrices.Mij2(i, i, hz, hz);
                                for (int j = 0; j < i; j++)
                                {
                                    int k = mat.ia[index[i]];
                                    while (mat.ja[k] < index[j])
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[j])
                                    {
                                        mat.al[k] += curitem.betta * Matrices.Mij2(i, j, hz, hz);
                                        mat.au[k] += curitem.betta * Matrices.Mij2(i, j, hz, hz);
                                    }
                                }
                            }

                        }
                        break;
                    case edges.Back:
                        foreach (var el in elems.FindAll((T) => T.ylownum + 1 == Ygrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 4), GetGlobalNum(el, 6), GetGlobalNum(el, 12), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 5), GetGlobalNum(el, 7), GetGlobalNum(el, 13), GetGlobalNum(el, 15) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double y = Ygrid[Ygrid.Count - 1];
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hx = x[1] - x[0];
                            for (int i = 0; i < 4; i++)
                            {
                                double locb = 0;
                                for (int j = 0; j < 4; j++)
                                {
                                    locb += curitem.betta * curitem.value(x[j % 2], y, z[j / 2]) * Matrices.Mij2(i, j, hz, hz);
                                }
                                mat.b[index[i]] += locb;
                                mat.di[index[i]] += curitem.betta * Matrices.Mij2(i, i, hz, hz);
                                for (int j = 0; j < i; j++)
                                {
                                    int k = mat.ia[index[i]];
                                    while (mat.ja[k] < index[j])
                                    {
                                        k++;
                                    }
                                    if (mat.ja[k] == index[j])
                                    {
                                        mat.al[k] += curitem.betta * Matrices.Mij2(i, j, hz, hz);
                                        mat.au[k] += curitem.betta * Matrices.Mij2(i, j, hz, hz);
                                    }
                                }
                            }
                        }
                        break;
                    default:
                        break;
                }
            }
        }
        public void AddBoundary2()
        {
            foreach (var item in BCs.FindAll((T) => T.bctype == 2))
            {
                var curitem = item as BC2;
                switch (curitem.edge)
                {
                    case edges.Top:
                        foreach (var el in elems.FindAll((T) => T.zlownum + 1 == Zgrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 8), GetGlobalNum(el, 10), GetGlobalNum(el, 12), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 9), GetGlobalNum(el, 11), GetGlobalNum(el, 13), GetGlobalNum(el, 15) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double z = Zgrid[Zgrid.Count - 1];
                            double hx = x[1] - x[0], hy = y[1] - y[0];
                            double[] locb = new double[4];
                            for (int i = 0; i < 4; i++)
                            {
                                for (int j = 0; j < 4; j++)
                                {
                                    locb[i] += curitem.value(x[j % 2], y[j / 2], z) * Matrices.Mij2(i, j, hx, hy);
                                }
                                mat.b[index[i]] += locb[i];
                            }
                        }
                        break;
                    case edges.Bottom:
                        foreach (var el in elems.FindAll((T) => T.zlownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 2), GetGlobalNum(el, 4), GetGlobalNum(el, 6) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 3), GetGlobalNum(el, 5), GetGlobalNum(el, 7) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double z = Zgrid[0];
                            double hx = x[1] - x[0], hy = y[1] - y[0];
                            double[] locb = new double[4];
                            for (int i = 0; i < 4; i++)
                            {
                                for (int j = 0; j < 4; j++)
                                {
                                    locb[i] += curitem.value(x[j % 2], y[j / 2], z) * Matrices.Mij2(i, j, hx, hy);
                                }
                                mat.b[index[i]] += locb[i];
                            }
                        }
                        break;
                    case edges.Left:
                        foreach (var el in elems.FindAll((T) => T.xlownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 4), GetGlobalNum(el, 8), GetGlobalNum(el, 12) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 5), GetGlobalNum(el, 9), GetGlobalNum(el, 13) };
                            }
                            double x = Xgrid[0];
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hy = y[1] - y[0];
                            double[] locb = new double[4];
                            for (int i = 0; i < 4; i++)
                            {
                                for (int j = 0; j < 4; j++)
                                {
                                    locb[i] += curitem.value(x, y[j % 2], z[j / 2]) * Matrices.Mij2(i, j, hz, hy);
                                }
                                mat.b[index[i]] += locb[i];
                            }
                        }
                        break;
                    case edges.Right:
                        foreach (var el in elems.FindAll((T) => T.xlownum + 1 == Xgrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 2), GetGlobalNum(el, 6), GetGlobalNum(el, 10), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 3), GetGlobalNum(el, 7), GetGlobalNum(el, 11), GetGlobalNum(el, 15) };
                            }
                            double x = Xgrid[Xgrid.Count - 1];
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hy = y[1] - y[0];
                            double[] locb = new double[4];
                            for (int i = 0; i < 4; i++)
                            {
                                for (int j = 0; j < 4; j++)
                                {
                                    locb[i] += curitem.value(x, y[j % 2], z[j / 2]) * Matrices.Mij2(i, j, hz, hy);
                                }
                                mat.b[index[i]] += locb[i];
                            }
                        }
                        break;
                    case edges.Front:
                        foreach (var el in elems.FindAll((T) => T.ylownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 2), GetGlobalNum(el, 8), GetGlobalNum(el, 10) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 3), GetGlobalNum(el, 9), GetGlobalNum(el, 11) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double y = Ygrid[0];
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hx = x[1] - x[0];
                            double[] locb = new double[4];
                            for (int i = 0; i < 4; i++)
                            {
                                for (int j = 0; j < 4; j++)
                                {
                                    locb[i] += curitem.value(x[j % 2], y, z[j / 2]) * Matrices.Mij2(i, j, hz, hx);
                                }
                                mat.b[index[i]] += locb[i];
                            }
                        }
                        break;
                    case edges.Back:
                        foreach (var el in elems.FindAll((T) => T.ylownum + 1 == Ygrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 4), GetGlobalNum(el, 6), GetGlobalNum(el, 12), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 5), GetGlobalNum(el, 7), GetGlobalNum(el, 13), GetGlobalNum(el, 15) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double y = Ygrid[Ygrid.Count - 1];
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            double hz = z[1] - z[0], hx = x[1] - x[0];
                            double[] locb = new double[4];
                            for (int i = 0; i < 4; i++)
                            {
                                for (int j = 0; j < 4; j++)
                                {
                                    locb[i] += curitem.value(x[j % 2], y, z[j / 2]) * Matrices.Mij2(i, j, hz, hx);
                                }
                                mat.b[index[i]] += locb[i];
                            }
                        }
                        break;
                    default:
                        break;
                }
            }
        }
        public void AddBoundary1()
        {
            foreach (var item in BCs.FindAll((T) => T.bctype == 1))
            {
                var curitem = item as BC1;
                switch (curitem.edge)
                {
                    case edges.Top:
                        foreach (var el in elems.FindAll((T) => T.zlownum + 1 == Zgrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 8), GetGlobalNum(el, 10), GetGlobalNum(el, 12), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 9), GetGlobalNum(el, 11), GetGlobalNum(el, 13), GetGlobalNum(el, 15) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double z = Zgrid[Zgrid.Count - 1];
                            for (int p = 0; p < 4; p++)
                            {
                                mat.di[index[p]] = 1;
                                mat.b[index[p]] = curitem.value(x[p % 2], y[p / 2], z);
                                for (int k = mat.ia[index[p]]; k < mat.ia[index[p] + 1]; k++)
                                {
                                    mat.b[mat.ja[k]] -= curitem.value(x[p % 2], y[p / 2], z) * mat.au[k];
                                    mat.al[k] = 0;
                                    mat.au[k] = 0;
                                }
                                for (int i = index[p] + 1; i < mat.n; i++)
                                {
                                    int low = mat.ia[i];
                                    int high = mat.ia[i + 1];
                                    bool flag = false;
                                    int mid = 0;
                                    while (low <= high)
                                    {
                                        mid = (low + high) / 2;
                                        int midVal = mat.ja[mid];
                                        if (midVal < index[p])
                                            low = mid + 1;
                                        else
                                        {
                                            if (midVal > index[p])
                                                high = mid - 1;
                                            else
                                            {
                                                flag = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (flag)
                                    {
                                        mat.b[i] -= mat.al[mid] * curitem.value(x[p % 2], y[p / 2], z);
                                        mat.au[mid] = 0;
                                        mat.al[mid] = 0;
                                    }
                                }
                            }
                        }
                        break;
                    case edges.Bottom:
                        foreach (var el in elems.FindAll((T) => T.zlownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 2), GetGlobalNum(el, 4), GetGlobalNum(el, 6) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 3), GetGlobalNum(el, 5), GetGlobalNum(el, 7) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double z = Zgrid[0];
                            for (int p = 0; p < 4; p++)
                            {
                                mat.di[index[p]] = 1;
                                mat.b[index[p]] = curitem.value(x[p % 2], y[p / 2], z);
                                for (int k = mat.ia[index[p]]; k < mat.ia[index[p] + 1]; k++)
                                {
                                    mat.b[mat.ja[k]] -= curitem.value(x[p % 2], y[p / 2], z) * mat.au[k];
                                    mat.au[k] = 0;
                                    mat.al[k] = 0;
                                }
                                for (int i = index[p] + 1; i < mat.n; i++)
                                {
                                    int low = mat.ia[i];
                                    int high = mat.ia[i + 1];
                                    bool flag = false;
                                    int mid = 0;
                                    while (low <= high)
                                    {
                                        mid = (low + high) / 2;
                                        int midVal = mat.ja[mid];
                                        if (midVal < index[p])
                                            low = mid + 1;
                                        else
                                        {
                                            if (midVal > index[p])
                                                high = mid - 1;
                                            else
                                            {
                                                flag = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (flag)
                                    {
                                        mat.b[i] -= mat.al[mid] * curitem.value(x[p % 2], y[p / 2], z);
                                        mat.au[mid] = 0;
                                        mat.al[mid] = 0;
                                    }
                                }
                            }
                        }
                        break;
                    case edges.Left:
                        foreach (var el in elems.FindAll((T) => T.xlownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 4), GetGlobalNum(el, 8), GetGlobalNum(el, 12) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 5), GetGlobalNum(el, 9), GetGlobalNum(el, 13) };
                            }
                            double x = Xgrid[0];
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            for (int p = 0; p < 4; p++)
                            {
                                mat.di[index[p]] = 1;
                                mat.b[index[p]] = curitem.value(x, y[p % 2], z[p / 2]);
                                for (int k = mat.ia[index[p]]; k < mat.ia[index[p] + 1]; k++)
                                {
                                    mat.b[mat.ja[k]] -= curitem.value(x, y[p % 2], z[p / 2]) * mat.au[k];
                                    mat.al[k] = 0;
                                    mat.au[k] = 0;
                                }
                                for (int i = index[p] + 1; i < mat.n; i++)
                                {
                                    int low = mat.ia[i];
                                    int high = mat.ia[i + 1];
                                    bool flag = false;
                                    int mid = 0;
                                    while (low <= high)
                                    {
                                        mid = (low + high) / 2;
                                        int midVal = mat.ja[mid];
                                        if (midVal < index[p])
                                            low = mid + 1;
                                        else
                                        {
                                            if (midVal > index[p])
                                                high = mid - 1;
                                            else
                                            {
                                                flag = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (flag)
                                    {
                                        mat.b[i] -= mat.al[mid] * curitem.value(x, y[p % 2], z[p / 2]);
                                        mat.au[mid] = 0;
                                        mat.al[mid] = 0;
                                    }
                                }
                            }
                        }
                        break;
                    case edges.Right:
                        foreach (var el in elems.FindAll((T) => T.xlownum + 1 == Xgrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 2), GetGlobalNum(el, 6), GetGlobalNum(el, 10), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 3), GetGlobalNum(el, 7), GetGlobalNum(el, 11), GetGlobalNum(el, 15) };
                            }
                            double x = Xgrid[Xgrid.Count - 1];
                            double[] y = { Ygrid[el.ylownum], Ygrid[el.ylownum + 1] };
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            for (int p = 0; p < 4; p++)
                            {
                                mat.di[index[p]] = 1;
                                mat.b[index[p]] = curitem.value(x, y[p % 2], z[p / 2]);
                                for (int k = mat.ia[index[p]]; k < mat.ia[index[p] + 1]; k++)
                                {
                                    mat.b[mat.ja[k]] -= curitem.value(x, y[p % 2], z[p / 2]) * mat.au[k];
                                    mat.al[k] = 0;
                                    mat.au[k] = 0;
                                }
                                for (int i = index[p] + 1; i < mat.n; i++)
                                {
                                    int low = mat.ia[i];
                                    int high = mat.ia[i + 1];
                                    bool flag = false;
                                    int mid = 0;
                                    while (low <= high)
                                    {
                                        mid = (low + high) / 2;
                                        int midVal = mat.ja[mid];
                                        if (midVal < index[p])
                                            low = mid + 1;
                                        else
                                        {
                                            if (midVal > index[p])
                                                high = mid - 1;
                                            else
                                            {
                                                flag = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (flag)
                                    {
                                        mat.b[i] -= mat.al[mid] * curitem.value(x, y[p % 2], z[p / 2]);
                                        mat.au[mid] = 0;
                                        mat.al[mid] = 0;
                                    }
                                }
                            }
                        }
                        break;
                    case edges.Front:
                        foreach (var el in elems.FindAll((T) => T.ylownum == 0))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 0), GetGlobalNum(el, 2), GetGlobalNum(el, 8), GetGlobalNum(el, 10) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 1), GetGlobalNum(el, 3), GetGlobalNum(el, 9), GetGlobalNum(el, 11) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double y = Ygrid[0];
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            for (int p = 0; p < 4; p++)
                            {
                                mat.di[index[p]] = 1;
                                mat.b[index[p]] = curitem.value(x[p % 2], y, z[p / 2]);
                                for (int k = mat.ia[index[p]]; k < mat.ia[index[p] + 1]; k++)
                                {
                                    mat.b[mat.ja[k]] -= curitem.value(x[p % 2], y, z[p / 2]) * mat.au[k];
                                    mat.al[k] = 0;
                                    mat.au[k] = 0;
                                }
                                for (int i = index[p] + 1; i < mat.n; i++)
                                {
                                    int low = mat.ia[i];
                                    int high = mat.ia[i + 1];
                                    bool flag = false;
                                    int mid = 0;
                                    while (low <= high)
                                    {
                                        mid = (low + high) / 2;
                                        int midVal = mat.ja[mid];
                                        if (midVal < index[p])
                                            low = mid + 1;
                                        else
                                        {
                                            if (midVal > index[p])
                                                high = mid - 1;
                                            else
                                            {
                                                flag = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (flag)
                                    {
                                        mat.b[i] -= mat.al[mid] * curitem.value(x[p % 2], y, z[p / 2]);
                                        mat.au[mid] = 0;
                                        mat.al[mid] = 0;
                                    }
                                }
                            }
                        }
                            break;
                    case edges.Back:
                        foreach (var el in elems.FindAll((T) => T.ylownum + 1 == Ygrid.Count - 1))
                        {
                            int[] index = null;
                            if (curitem.cossintype == 0)//cos
                            {
                                index = new int[] { GetGlobalNum(el, 4), GetGlobalNum(el, 6), GetGlobalNum(el, 12), GetGlobalNum(el, 14) };

                            }
                            else//sin
                            {
                                index = new int[] { GetGlobalNum(el, 5), GetGlobalNum(el, 7), GetGlobalNum(el, 13), GetGlobalNum(el, 15) };
                            }
                            double[] x = { Xgrid[el.xlownum], Xgrid[el.xlownum + 1] };
                            double y = Ygrid[Ygrid.Count - 1];
                            double[] z = { Zgrid[el.zlownum], Zgrid[el.zlownum + 1] };
                            for (int p = 0; p < 4; p++)
                            {
                                mat.di[index[p]] = 1;
                                mat.b[index[p]] = curitem.value(x[p % 2], y, z[p / 2]);
                                for (int k = mat.ia[index[p]]; k < mat.ia[index[p] + 1]; k++)
                                {
                                    mat.b[mat.ja[k]] -= curitem.value(x[p % 2], y, z[p / 2]) * mat.au[k];
                                    mat.al[k] = 0;
                                    mat.au[k] = 0;
                                }
                                for (int i = index[p] + 1; i < mat.n; i++)
                                {
                                    int low = mat.ia[i];
                                    int high = mat.ia[i + 1];
                                    bool flag = false;
                                    int mid = 0;
                                    while (low <= high)
                                    {
                                        mid = (low + high) / 2;
                                        int midVal = mat.ja[mid];
                                        if (midVal < index[p])
                                            low = mid + 1;
                                        else
                                        {
                                            if (midVal > index[p])
                                                high = mid - 1;
                                            else
                                            {
                                                flag = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (flag)
                                    {
                                        mat.b[i] -= mat.al[mid] * curitem.value(x[p % 2], y, z[p / 2]);
                                        mat.au[mid] = 0;
                                        mat.al[mid] = 0;
                                    }
                                }
                            }
                        }
                        break;
                    default:
                        break;
                }
            }
        }
        public void SolveLU()
        {
            GenerateProfile();
            foreach (var item in elems)
            {
                Addlocal(item);
            }
            AddBoundary2();
            AddBoundary3();
            AddBoundary1();
            matprof = mat.ToProfile();
            q = matprof.Solve();
        }
        public void SolveLOS(double eps, int maxiter)
        {
            GenerateProfile();
            foreach (var item in elems)
            {
                Addlocal(item);
            }
            AddBoundary2();
            AddBoundary3();
            AddBoundary1();
            List<double> x0 = new();
            for (int i = 0; i < mat.n; i++)
            {
                x0.Add(0);
            }
            q = mat.LOS(x0, maxiter, eps);
        }
        public void SolveLOSPrecond(double eps, int maxiter)
        {
            GenerateProfile();
            foreach (var item in elems)
            {
                Addlocal(item);
            }
            AddBoundary2();
            AddBoundary3();
            AddBoundary1();
            List<double> x0 = new();
            for (int i = 0; i < mat.n; i++)
            {
                x0.Add(0);
            }
            mat.LU();
            q = mat.LoS_precond(x0, eps, maxiter);
        }
        public void SolveGMRES(double eps, int maxiter, int depth)
        {
            GenerateProfile();
            foreach (var item in elems)
            {
                Addlocal(item);
            }
            AddBoundary2();
            AddBoundary3();
            AddBoundary1();
            List<double> x0 = new();
            for (int i = 0; i < mat.n; i++)
            {
                x0.Add(0);
            }
            mat.LU();
            q = mat.GMRES(x0, eps, maxiter, depth);
        }
        public (double, double) GetSollution(double x, double y, double z)
        {
            bool flag = true;
            int n = elems.Count;
            int i = 0;
            for (i = 0; i < n && flag; i++)
            {
                if (x >= Xgrid[elems[i].xlownum] && x <= Xgrid[elems[i].xlownum + 1] && y >= Ygrid[elems[i].ylownum] && y <= Ygrid[elems[i].ylownum + 1] && z >= Zgrid[elems[i].zlownum] && z <= Zgrid[elems[i].zlownum + 1])
                {
                    flag = false;
                    i--;
                }
            }
            if (flag)
                throw new Exception();
            double xloc = (x - Xgrid[elems[i].xlownum]) / (Xgrid[elems[i].xlownum + 1] - Xgrid[elems[i].xlownum]);
            double yloc = (y - Ygrid[elems[i].ylownum]) / (Ygrid[elems[i].ylownum + 1] - Ygrid[elems[i].ylownum]);
            double zloc = (z - Zgrid[elems[i].zlownum]) / (Zgrid[elems[i].zlownum + 1] - Zgrid[elems[i].zlownum]);
            double[] psix = { 1 - xloc, xloc };
            double[] psiy = { 1 - yloc, yloc };
            double[] psiz = { 1 - zloc, zloc };
            (double, double) res = (0, 0);
            for (int p = 0; p < 8; p++)
            {
                res.Item1 += q[GetGlobalNum(elems[i], 2 * p)] * psix[p % 2] * psiy[(p / 2) % 2] * psiz[p / 4];
                res.Item2 += q[GetGlobalNum(elems[i], 2 * p + 1)] * psix[p % 2] * psiy[(p / 2) % 2] * psiz[p / 4];
            }
            return res;
        }
    }
    public class Elem
    {
        public Elem(int xlownum, int ylownum, int zlownum)
        {
            this.xlownum = xlownum;
            this.ylownum = ylownum;
            this.zlownum = zlownum;
        }

        public int xlownum { get; init; }
        public int ylownum { get; init; }
        public int zlownum { get; init; }
    }
    public abstract class BC
    {
        public edges edge { get; init; }
        public int cossintype { get; init; }
        public int bctype { get; init; }

        public Func<double, double, double, double> value;

        public BC(int bctype, edges edge, Func<double, double, double, double> value, int cossintype)
        {
            this.bctype = bctype;
            this.edge = edge;
            this.value = value;
            this.cossintype = cossintype;
        }
    }
    public class BC1 : BC
    {
        public BC1(edges edge, Func<double, double, double, double> value, int cossintype) : base(1, edge, value, cossintype) { }
    }
    public class BC2 : BC
    {
        public BC2(edges edge, Func<double, double, double, double> value, int cossintype) : base(2, edge, value, cossintype) { }
    }
    public class BC3 : BC
    {
        public double betta { get; init; }
        public BC3(edges edge, Func<double, double, double, double> value, double betta, int cossintype) : base(3, edge, value, cossintype) { this.betta = betta; }
    }
    public static class Matrices
    {
        static double[][] Gmatr = new double[][]
        {
            new double[]
            {
                1,-1
            },

            new double[]
            {
                -1,1
            }
        };
        static double[][] Mmatr = new double[][]
       {
            new double[]
            {
               2.0/6,1.0/6
            },
            new double[]
            {
                1.0/6,2.0/6
            }
       };
        public static double Gij(int i, int j, double hx, double hy, double hz)
        {
            return hy * hz / hx * Gmatr[i % 2][j % 2] * Mmatr[(i / 2) % 2][(j / 2) % 2] * Mmatr[i / 4][j / 4] + hx * hz / hy * Mmatr[i % 2][j % 2] * Gmatr[(i / 2) % 2][(j / 2) % 2] * Mmatr[i / 4][j / 4] + hx * hy / hz * Mmatr[i % 2][j % 2] * Mmatr[(i / 2) % 2][(j / 2) % 2] * Gmatr[i / 4][j / 4];
        }
        public static double Mij(int i, int j, double hx, double hy, double hz) =>
            hx * hy * hz * Mmatr[i % 2][j % 2] * Mmatr[(i / 2) % 2][(j / 2) % 2] * Mmatr[i / 4][j / 4];
        public static double Mij2(int i, int j, double hx, double hy) =>
            hx * hy * Mmatr[i % 2][j % 2] * Mmatr[i / 2][j / 2];
    }
    public class MatrixSparce
    {
        public List<double> al;
        public List<double> au;
        public List<double> di;
        public List<int> ia;
        public List<int> ja;
        public List<double> b;
        public List<double> di_LU = new();
        public List<double> au_LU = new();
        public List<double> al_LU = new();
        public int n;
        public MatrixSparce() { }
        public List<double> MSG(List<double> x0, int maxiter, double eps)
        {
            double bnorm = Math.Sqrt(DotProduct(b, b));
            List<double> p = new List<double>(n);
            List<double> q = new List<double>(n);
            for (int i = 0; i < n; i++)
            {
                p.Add(0);
            }
            var r = MatrixMult(x0);
            for (int i = 0; i < n; i++)
            {
                r[i] = b[i] - r[i];
                p[i] = r[i];
            }
            int k = 0;
            double alpha, betta, rnorm = Math.Sqrt(DotProduct(r, r));
            while (k < maxiter && rnorm / bnorm > eps)
            {
                q = MatrixMult(p);
                alpha = DotProduct(r, r) / DotProduct(q, p);
                betta = 1 / DotProduct(r, r);
                for (int i = 0; i < n; i++)
                {
                    x0[i] += alpha * p[i];
                    r[i] -= alpha * q[i];
                }
                rnorm = Math.Sqrt(DotProduct(r, r));
                betta *= DotProduct(r, r);
                for (int i = 0; i < n; i++)
                {
                    p[i] = r[i] + betta * p[i];
                }
                k++;
            }
            return x0;
        }
        public void LU()
        {

            foreach (var item in di)
            {
                di_LU.Add(item);
            }
            foreach (var item in au)
            {
                au_LU.Add(item);
            }
            foreach (var item in al)
            {
                al_LU.Add(item);
            }
            for (int i = 0; i < n; i++)
            {
                double sumdi = 0.0;

                int i0 = ia[i];
                int i1 = ia[i + 1];


                for (int k = i0; k < i1; k++)
                {
                    int j = ja[k];
                    int j0 = ia[j];

                    int j1 = ia[j + 1];


                    int ik = i0;
                    int kj = j0;

                    double suml = 0.0;
                    double sumu = 0.0;

                    while (ik < k)
                    {

                        if (ja[ik] == ja[kj])
                        {

                            suml += al_LU[ik] * au_LU[kj];
                            sumu += au_LU[ik] * al_LU[kj];
                            ik++;
                            kj++;
                        }

                        else
                        {
                            if (ja[ik] > ja[kj])
                            {
                                kj++;
                            }
                            else
                            {
                                ik++;
                            }
                        }
                    }

                    al_LU[k] = al_LU[k] - suml;
                    au_LU[k] = (au_LU[k] - sumu) / di_LU[j];
                    sumdi += al_LU[k] * au_LU[k];
                }

                di_LU[i] = di_LU[i] - sumdi;
            }
        }
        public List<double> LoS_precond(List<double> x0, double eps, int maxiter)
        {
            int k = 1;
            List<double> buf = MatrixMult(x0);
            double bnorm = 0;
            for (int i = 0; i < n; i++)
            {
                buf[i] = b[i] - buf[i];
            }
            double rnorm = Math.Sqrt(DotProduct(buf, buf));
            List<double> r = LUDirect(buf);
            bnorm = Math.Sqrt(DotProduct(b, b));
            List<double> z = LUReverse(r);
            buf = MatrixMult(z);
            List<double> p = LUDirect(buf);
            double resid = 1;
            while (resid > eps && rnorm / bnorm > eps * eps && k < maxiter)
            {
                double pp = DotProduct(p, p);
                double pr = DotProduct(p, r);
                double alpha = pr / pp;
                for (int i = 0; i < n; i++)
                {
                    x0[i] += alpha * z[i];
                    r[i] -= alpha * p[i];
                }
                rnorm = Math.Sqrt(DotProduct(r, r));
                List<double> Ur = LUReverse(r);
                buf = MatrixMult(Ur);
                buf = LUDirect(buf);
                double betta = -(DotProduct(p, buf) / pp);
                for (int i = 0; i < n; i++)
                {
                    z[i] = Ur[i] + betta * z[i];
                    p[i] = buf[i] + betta * p[i];
                }
                double test1 = 0;
                double test2 = 0;
                var asd = MatrixMult(x0);
                for (int i = 0; i < n; i++)
                {
                    test1 += (asd[i] - b[i]) * (asd[i] - b[i]);
                    test2 += b[i] * b[i];
                }
                resid = Math.Sqrt(test1 / test2);
                k++;
            }
            Console.WriteLine($"{k} {rnorm / bnorm} {resid}");
            return x0;
        }
        List<double> LUDirect(List<double> rpart)
        {
            List<double> res = new();
            for (int i = 0; i < n; i++)
            {
                res.Add(rpart[i]);
            }

            for (int i = 0; i < n; i++)
            {
                double sum = 0.0;
                for (int j = ia[i]; j < ia[i + 1]; j++)
                    sum += al_LU[j] * res[ja[j]];
                res[i] -= sum;
                res[i] /= di_LU[i];
            }
            return res;
        }
        List<double> LUReverse(List<double> rpart)
        {
            List<double> res = new();
            for (int i = 0; i < n; i++)
            {
                res.Add(rpart[i]);
            }

            for (int i = n - 1; i >= 0; i--)
            {
                for (int j = ia[i]; j < ia[i + 1]; j++)
                    res[ja[j]] -= au_LU[j] * res[i];
            }
            return res;
        }
        List<double> MultU(List<double> x0)
        {
            List<double> res = new();
            for (int i = 0; i < n; i++)
            {
                res.Add(x0[i]);
                int k0 = ia[i], k1 = ia[i + 1];
                for (int k = k0; k < k1; k++)
                {
                    res[i] += au_LU[k] * x0[ja[k]];
                }
            }
            return res;
        }
        public List<double> LOS(List<double> x0, int maxiter, double eps)
        {
            double bnorm = Math.Sqrt(DotProduct(b, b));

            List<double> Ar = new();
            List<double> r = MatrixMult(x0);
            List<double> z = new();
            for (int i = 0; i < n; i++)
            {
                r[i] = b[i] - r[i];
                z.Add(r[i]);
            }
            List<double> p = MatrixMult(z);
            int k = 0;
            double alpha, betta, rnorm = Math.Sqrt(DotProduct(r, r));
            while (k < maxiter && rnorm / bnorm > eps)
            {
                alpha = DotProduct(p, r) / DotProduct(p, p);
                for (int i = 0; i < n; i++)
                {
                    x0[i] += alpha * z[i];
                    r[i] -= alpha * p[i];
                }
                Ar = MatrixMult(r);
                betta = -DotProduct(p, Ar) / DotProduct(p, p);
                rnorm = Math.Sqrt(DotProduct(r, r));
                for (int i = 0; i < n; i++)
                {
                    z[i] = r[i] + betta * z[i];
                    p[i] = Ar[i] + betta * p[i];
                }
                k++;
            }
            Console.WriteLine($"{rnorm / bnorm} {k}");
            return x0;
        }
        public List<double> GMRES(List<double> x0, double eps, int maxiter, int depth)
        {
            //Q верхнетреугольная обратный 
            //S нижнетреугольная прямой
            int curdepth = depth;
            double bnorm = Math.Sqrt(DotProduct(b, b));
            List<double> x = MultU(x0);//меняем на x с волной

            var r = MatrixMult(x0);
            for (int i = 0; i < n; i++)
            {
                r[i] = b[i] - r[i];
            }
            r = LUDirect(r);//считаем невязку с волной

            List<List<double>> V = new();//транспонированная
            double[][] H = new double[depth + 1][];
            for (int i = 0; i < depth + 1; i++)
            {
                H[i] = new double[depth];
            }
            double rnorm = Math.Sqrt(DotProduct(r, r));
            int iter = 0;
            while (rnorm/bnorm>eps&&iter<maxiter)
            {
                curdepth = depth;
                V.Add(new());
                for (int i = 0; i < n; i++)
                {
                    V[0].Add(r[i] / rnorm);// добавляем первый вектор в матрицу V
                }
                for (int i = 0; i < depth; i++)// заполняется матрица H
                {
                    List<double> w = LUReverse(V[i]);
                    w = MatrixMult(w);
                    w = LUDirect(w);
                    for (int j = 0; j <= i; j++)
                    {
                        H[j][i] = DotProduct(V[j], w);
                    }
                    for (int j = 0; j <= i; j++)
                    {
                        for (int k = 0; k < n; k++)
                        {
                            w[k] -= H[j][i] * V[j][k];
                        }
                    }

                    H[i + 1][i] = Math.Sqrt(DotProduct(w, w));
                    if (H[i + 1][i] == 0)// если новый вектор нулевой заканчиваем заполнение матриц
                    {
                        curdepth = i + 1;
                        break;
                    }
                    else
                    {
                        V.Add(new());
                        for (int k = 0; k < n; k++)
                        {
                            V[i + 1].Add(w[k] / H[i + 1][i]);//если ненулевой добавляем в матрицу V
                        }
                    }
                }
                List<double> d = new();//дальше будет минимизация d-Hz
                d.Add(rnorm);
                for (int i = 1; i < curdepth + 1; i++)
                {
                    d.Add(0);
                }
                for (int i = 0; i < curdepth; i++)//умножаем правую часть и матрицу на матрицу поворота чтобы убрать диагональ под главной
                {
                    double norm = Math.Sqrt(H[i][i] * H[i][i] + H[i + 1][i] * H[i + 1][i]);
                    double c = H[i][i] / norm;
                    double s = H[i + 1][i] / norm;
                    for (int k = i; k < curdepth; k++)
                    {
                        double ii = c * H[i][k] + s * H[i + 1][k];
                        double i1i = -s * H[i][k] + c * H[i + 1][k];
                        H[i][k] = ii;
                        H[i + 1][k] = i1i;
                    }
                    double d1 = d[i] * c + d[i + 1] * s;
                    double d2 = -s * d[i] + c * d[i + 1];
                    d[i] = d1;
                    d[i + 1] = d2;
                }

                for (int i = curdepth - 1; i >= 0; i--)//обратный гаус для матрицы H верхнетреугольной
                {
                    double summ = 0;
                    for (int j = curdepth - 1; j > i; j--)
                    {
                        summ += d[j] * H[i][j];
                    }
                    d[i] = (d[i] - summ) / H[i][i];
                }

                for (int i = 0; i < n; i++)//делаем добавку в x
                {
                    for (int j = 0; j < curdepth; j++)
                    {
                        x[i] += V[j][i] * d[j];
                    }
                }

                r = LUReverse(x);
                r = MatrixMult(r);
                for (int i = 0; i < n; i++)
                {
                    r[i] = b[i] - r[i];
                }
                r = LUDirect(r);
                rnorm = Math.Sqrt(DotProduct(r, r));
                iter++;
                V.Clear();
            }
            Console.WriteLine($"{iter} {rnorm / bnorm}");
            return LUReverse(x);
        }
        public void PrintDense()
        {
            for (int i = 0; i < n; i++)
            {
                int k = ia[i];
                for (int j = 0; j < i; j++)
                {
                    while (ja[k] < j && k < ia[i + 1] - 1)
                        k++;
                    if (ja[k] == j)
                        Console.Write($"{al[k]:f2} ");
                    else
                        Console.Write($"{0:f2} ");
                }
                Console.Write($"{di[i]:f2} ");
                for (int j = i + 1; j < n; j++)
                {
                    k = ia[j];
                    while (ja[k] < i && k < ia[i + 1] - 1)
                        k++;
                    if (ja[k] == i)
                        Console.Write($"{au[k]:f2} ");
                    else
                        Console.Write($"{0:f2} ");
                }

                Console.WriteLine();
            }
        }
        public MatrixProfile ToProfile()
        {
            MatrixProfile res = new();
            res.n = n;
            res.di = di;
            res.b = b;
            res.ia = new();
            res.al = new();
            res.au = new();
            res.ia.Add(0);
            res.ia.Add(0);
            for (int i = 1; i < n; i++)
            {
                if (ia[i + 1] - ia[i] > 0)
                {
                    int count = 0;
                    int k = ia[i];
                    for (int p = ja[ia[i]]; p < i; p++)
                    {
                        count++;
                        while (ja[k] < p)
                            k++;
                        if (ja[k] == p)
                        {
                            res.al.Add(al[k]);
                            res.au.Add(au[k]);
                        }
                        else
                        {
                            res.al.Add(0);
                            res.au.Add(0);
                        }
                    }
                    res.ia.Add(res.ia[res.ia.Count - 1] + count);
                }
                else
                {
                    res.ia.Add(res.ia[i - 1]);
                }
            }
            return res;
        }
        double DotProduct(List<double> vec1, List<double> vec2)
        {
            if (vec1.Count != vec2.Count)
                throw new Exception();
            double res = 0;
            int m = vec1.Count;
            for (int i = 0; i < m; i++)
            {
                res += vec1[i] * vec2[i];
            }
            return res;
        }
        List<double> MatrixMult(List<double> x)
        {
            if (x.Count != n)
                throw new Exception();
            List<double> res = new List<double>(x.Count);
            for (int i = 0; i < n; i++)
            {
                res.Add(0);
            }
            for (int i = 0; i < n; i++)
            {
                res[i] = x[i] * di[i];
                for (int k = ia[i]; k < ia[i + 1]; k++)
                {
                    int j = ja[k];
                    res[i] += al[k] * x[j];
                    res[j] += au[k] * x[i];
                }
            }
            return res;
        }
    }
    public class MatrixProfile
    {
        public List<double> al;
        public List<double> au;
        public List<double> di;
        public List<int> ia;
        public List<double> b;
        public int n;
        public MatrixProfile() { }
        void LU()
        {
            for (int i = 0; i < n; i++)
            {
                double sumrow = 0;
                int j0 = i - (ia[i + 1] - ia[i]);
                int dind = j0;
                for (int ii = ia[i]; ii < ia[i + 1]; ii++, dind++)
                {
                    int j = j0 + ii - ia[i];
                    int j0j = j - (ia[j + 1] - ia[j]);

                    int kbeg = j0 > j0j ? j0 : j0j;
                    int kend = i > j ? j : i;
                    double summu = 0;
                    double summl = 0;
                    int indexi = ia[i] + kbeg - j0;
                    int indexj = ia[j] + kbeg - j0j;
                    for (int k = 0; k < kend - kbeg; k++, indexi++, indexj++)
                    {
                        summu += al[indexj] * au[indexi];
                        summl += al[indexi] * au[indexj];
                    }
                    al[ii] -= summl;
                    au[ii] = (au[ii] - summu) / di[j];
                    sumrow += au[ii] * al[ii];
                }
                di[i] -= sumrow;
            }
        }
        void GaussU()
        {
            for (int i = n - 1; i >= 0; i--)
            {
                double cursum = 0;
                for (int j = n - 1; j > i; j--)
                {
                    if (ia[j + 1] - ia[j] >= j - i)
                    {
                        cursum += au[ia[j + 1] - j + i] * b[j];
                    }
                }
                b[i] -= cursum;
            }
        }
        void GaussL()
        {
            for (int i = 0; i < n; i++)
            {
                double cursum = 0;
                int j = i - (ia[i + 1] - ia[i]);
                for (int k = ia[i]; k < ia[i + 1]; k++, j++)
                {
                    cursum += al[k] * b[j];
                }
                b[i] = (b[i] - cursum) / di[i];
            }
        }
        public List<double> Solve()
        {
            LU();
            GaussL();
            GaussU();
            return b;
        }
    }

}
