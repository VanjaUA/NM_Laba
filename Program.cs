using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NM_Laba1
{
    class MainClass
    {
        public static void Main()
        {
            //GaussMethod.SolveGauss();
            //ChordMethod.ChordMethodSolver();
            //NewtonMethod.Run();
            //TrapezoidalMethod.Run();
            //RCL.Run();
            Laba6.Run();
        }
    }

    class GaussMethod
    {
        public static void SolveGauss()
        {
            double[,] matrix = {
            {8.3, 2.7, 4.1, 1.9, -10.15},
            {3.92, 8.45, 7.7, 2.46, 12.21},
            {3.77, 7.29, 8.04, 2.28, 15.4},
            {2.21, 3.57, 1.69, 6.69, -8.35},
            };

            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1) - 1;

            double[] solution = new double[cols];

            // Прямий хід
            for (int i = 0; i < rows - 1; i++)
            {
                // Вибір головного елемента по всій матриці
                int pivotRow = i;
                int pivotCol = i;
                double maxElement = Math.Abs(matrix[i, i]);

                for (int row = i; row < rows; row++)
                {
                    for (int col = i; col < cols; col++)
                    {
                        double currentElement = Math.Abs(matrix[row, col]);
                        if (currentElement > maxElement)
                        {
                            maxElement = currentElement;
                            pivotRow = row;
                            pivotCol = col;
                        }
                    }
                }

                // Обмін рядків
                for (int k = 0; k <= cols; k++)
                {
                    double temp = matrix[i, k];
                    matrix[i, k] = matrix[pivotRow, k];
                    matrix[pivotRow, k] = temp;
                }

                // Елімінація
                for (int j = i + 1; j < rows; j++)
                {
                    double factor = matrix[j, i] / matrix[i, i];
                    for (int k = i; k <= cols; k++)
                    {
                        matrix[j, k] -= factor * matrix[i, k];
                    }
                }
            }

            // Зворотній хід
            for (int i = rows - 1; i >= 0; i--)
            {
                solution[i] = matrix[i, cols] / matrix[i, i];
                for (int j = i - 1; j >= 0; j--)
                {
                    matrix[j, cols] -= matrix[j, i] * solution[i];
                }
            }

            Console.WriteLine("Результати:");
            for (int i = 0; i < solution.Length; i++)
            {
                Console.WriteLine($"x{i + 1} = {solution[i]}");
            }
            Console.ReadLine();
        }
    }

    class ChordMethod
    {
        public static void ChordMethodSolver()
        {
            double a = -1.5;
            double b = 0;
            double epsilon = 0.001;
            int maxIterations = 100;


            double Equation(double param)
            {
                return Math.Cos(param) - (1 / (param + 2));
            }



            double x = a;
            double result;
            for (int i = 0; i < maxIterations; i++)
            {
                double xOld = x;
                double fa = Equation(a);
                double fb = Equation(b);

                x = a - fa * ((b - a) / (fb - fa));

                if (Equation(x) * fa > 0)
                {
                    a = x;
                }
                else
                {
                    b = x;
                }

                if (Math.Abs((x - xOld) / x) * 100 < epsilon)
                {
                    result = x;
                    break;
                }
            }
            result = x;
            Console.WriteLine($"Корінь рівняння: {result}");
            Console.ReadLine();
        }
    }

    class NewtonMethod
    {
        // Делегат для представлення системи нелінійних рівнянь
        public delegate double[] SystemEquationsDelegate(double[] xy);

        // Похідні частинні за x та y
        public delegate double[,] JacobianMatrixDelegate(double[] xy);

        // Метод обертання матриці методом Гауса з вибором головного елемента по рядку
        private static void GaussianEliminationWithPartialPivoting(double[,] augmentedMatrix)
        {
            int n = augmentedMatrix.GetLength(0);

            for (int i = 0; i < n; i++)
            {
                // Знаходження максимального елементу в стовпці, що стоїть на діагоналі
                double maxElement = Math.Abs(augmentedMatrix[i, i]);
                int maxRow = i;
                for (int k = i + 1; k < n; k++)
                {
                    if (Math.Abs(augmentedMatrix[k, i]) > maxElement)
                    {
                        maxElement = Math.Abs(augmentedMatrix[k, i]);
                        maxRow = k;
                    }
                }

                // Обмін рядками для забезпечення ненульового головного елемента
                if (maxRow != i)
                {
                    for (int k = i; k < n + 1; k++)
                    {
                        double temp = augmentedMatrix[i, k];
                        augmentedMatrix[i, k] = augmentedMatrix[maxRow, k];
                        augmentedMatrix[maxRow, k] = temp;
                    }
                }

                // Обчислення коефіцієнтів для обнулення елементів під діагоналлю
                for (int k = i + 1; k < n; k++)
                {
                    double factor = augmentedMatrix[k, i] / augmentedMatrix[i, i];
                    for (int j = i; j < n + 1; j++)
                    {
                        augmentedMatrix[k, j] -= factor * augmentedMatrix[i, j];
                    }
                }
            }

            // Обчислення розв'язку системи за зворотнім ходом
            for (int i = n - 1; i >= 0; i--)
            {
                augmentedMatrix[i, n] /= augmentedMatrix[i, i];

                for (int k = i - 1; k >= 0; k--)
                {
                    augmentedMatrix[k, n] -= augmentedMatrix[k, i] * augmentedMatrix[i, n];
                }
            }
        }

        // Спрощений метод Ньютона
        public static double[] NewtonSimplified(SystemEquationsDelegate systemEquations, JacobianMatrixDelegate jacobianMatrix, double[] initialGuess, double epsilon, int maxIterations)
        {
            int n = initialGuess.Length;
            double[] xy = new double[n];

            Array.Copy(initialGuess, xy, n);

            for (int iteration = 0; iteration < maxIterations; iteration++)
            {
                double[] f = systemEquations(xy);
                double[,] jacobian = jacobianMatrix(xy);

                double[,] augmentedMatrix = new double[n, n + 1];

                // Створення розширеної матриці [J | -F]
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        augmentedMatrix[i, j] = jacobian[i, j];
                    }
                    augmentedMatrix[i, n] = -f[i];
                }

                // Виклик методу Гауса для розв'язання системи лінійних рівнянь
                GaussianEliminationWithPartialPivoting(augmentedMatrix);

                // Оновлення значень xy
                for (int i = 0; i < n; i++)
                {
                    xy[i] += augmentedMatrix[i, n];
                }

                // Перевірка на збіжність
                if (MaxNorm(f) < epsilon)
                {
                    return xy;
                }
            }

            return xy;
        }

        // Обчислення максимальної норми вектора
        private static double MaxNorm(double[] vector)
        {
            double maxNorm = Math.Abs(vector[0]);
            for (int i = 1; i < vector.Length; i++)
            {
                double absValue = Math.Abs(vector[i]);
                if (absValue > maxNorm)
                {
                    maxNorm = absValue;
                }
            }
            return maxNorm;
        }

        public static void Run()
        {
            //система рівнянь { 4*x^2 + y^2 - 4 = 0, x - y^2 + 1 = 0 }
            SystemEquationsDelegate system = xy => new double[]
            {
            4 * xy[0] * xy[0] + xy[1] * xy[1] - 4,
            xy[0] - xy[1] * xy[1] + 1
            };

            // Похідні частинні за x та y
            JacobianMatrixDelegate jacobian = xy => new double[,]
            {
            { 8 * xy[0], 2 * xy[1] },
            { 1, -2 * xy[1] }
            };

            double[] initialGuess = { 2, 2 }; // Початкові наближення для x і y
            double epsilon = 0.00001; //параметр точності
            int maxIterations = 100;

            double[] result = NewtonSimplified(system, jacobian, initialGuess, epsilon, maxIterations);

            Console.WriteLine($"Корені системи рівнянь: x1 = {result[0]}, х2 = {result[1]}");
            Console.ReadLine();
        }
    }

    class TrapezoidalMethod
    {
        static double Function(double x)
        {
            return Math.Pow(Math.Sin(x), 2);
        }

        static double TrapezoidalRule(double a, double b, int n)
        {
            double h = (b - a) / n; // Крок

            double result = (Function(a) + Function(b)) / 2.0;

            for (int i = 1; i < n; i++)
            {
                double x_i = a + i * h;
                result += Function(x_i);
            }

            result *= h;

            return result;
        }

        public static void Run()
        {
            double a = 0.0; // Нижня межа інтегрування
            double b = Math.PI; // Верхня межа інтегрування
            int n = 30;   

            double result = TrapezoidalRule(a, b, n);

            Console.WriteLine($"Значення інтегралу за допомогою методу трапецій: {result}");
            Console.ReadLine();
        }
    }

    class RCL 
    {
        static double t0 = 0;
        static double tn = 0.2;
        static double h = 1e-5;
        static double[] x0 = { 0, 0, 0 };
        public static void Run()
        {
            var t_values = GenerateTimeValues(t0, tn, h);
            var x_values = EulerMethodSystem(F, x0, t_values);

            var u2_values = new double[x_values.Length];
            var u1_values = new double[t_values.Length];
            for (int i = 0; i < x_values.Length; i++)
            {
                u2_values[i] = 4 * x_values[i][2];
                u1_values[i] = 100 * Math.Sin(2 * Math.PI * 50 * t_values[i]);
            }

            Console.WriteLine("Time (t)    u1        u2");
            for (int i = 0; i < t_values.Length; i++)
            {
                Console.WriteLine($"{t_values[i],-10:F6} {u1_values[i],-10:F6} {u2_values[i],-10:F6}");
            }
            Console.ReadLine();
        }

        static double[][] EulerMethodSystem(Func<double, double[], double[]> F, double[] x0, double[] t_values)
        {
            int numSteps = t_values.Length - 1;
            double[][] x_values = new double[t_values.Length][];
            x_values[0] = x0;

            for (int i = 0; i < numSteps; i++)
            {
                double[] x = x_values[i];
                double t = t_values[i];
                double[] x_next = new double[x.Length];
                double[] dxdt = F(t, x);

                for (int j = 0; j < x.Length; j++)
                {
                    x_next[j] = x[j] + h * dxdt[j];
                }

                x_values[i + 1] = x_next;
            }

            return x_values;
        }

        static double[] GenerateTimeValues(double t0, double tn, double h)
        {
            int numSteps = (int)((tn - t0) / h) + 1;
            double[] t_values = new double[numSteps];
            for (int i = 0; i < numSteps; i++)
            {
                t_values[i] = t0 + i * h;
            }
            return t_values;
        }

        static double[] F(double t, double[] x)
        {
            double dxdt = Math.Sin(2 * Math.PI * 50 * t) - x[0] - x[1] / 5 * 300e-6;
            double dydt = (Math.Sin(2 * Math.PI * 50 * t) - x[0] - x[1] / 5 - x[1] / 4 - (x[1] - x[2] * 4) * 7 / (4 + 7) * 4) * 1 / 150e-6;
            double dzdt = ((x[1] - x[2] * 4) / (4 * 7) * 0.02) * 7;

            return new double[] { dxdt, dydt, dzdt };
        }

    }



    class Laba6
    {
        static double R1 = 50;  // Ом
        static double R2 = 60;  // Ом
        static double R3 = 30;  // Ом
        static double C1 = 1.54;  // мФ
        static double C2 = 1;  // мФ
        static double L_min = 4.7;  // Гн
        static double L_max = 47;  // Гн
        static double i_min = 1;  // A
        static double i_max = 2;  // A

        // Алгоритм методу Рунге-Кутта для системи
        static void RungeKuttaSystem(double x0, double[] y0, double xEnd, double h)
        {
            double x = x0;
            double y1 = y0[0];
            double y2 = y0[1];
            double y3 = y0[2];

            while (x < xEnd)
            {
                Console.WriteLine($"час = {x}, " +
                                  $"напруга на конденсаторі с1 = {y1}, " +
                                  $"струм в індуктивності = {y2}, " +
                                  $"напруга на конденсаторі с2 = {y3}," +
                                  $"напруга u2 ={50 * y3}," +
                                  $"індуктивність l2={L2(y2)}," +
                                  $"напруга u1={u1(x)}");

                double k1 = h * F1(x, y1, y2, y3);
                double l1 = h * F2(x, y1, y2, y3);
                double m1 = h * F3(x, y1, y2, y3);

                double k2 = h * F1(x + h / 2, y1 + k1 / 2, y2 + l1 / 2, y3 + m1 / 2);
                double l2 = h * F2(x + h / 2, y1 + k1 / 2, y2 + l1 / 2, y3 + m1 / 2);
                double m2 = h * F3(x + h / 2, y1 + k1 / 2, y2 + l1 / 2, y3 + m1 / 2);

                double k3 = h * F1(x + h / 2, y1 + 2 * k2 - k1, y2 + 2 * l2 - l1, y3 + 2 * m2 - m1);
                double l3 = h * F2(x + h / 2, y1 + 2 * k2 - k1, y2 + 2 * l2 - l1, y3 + 2 * m2 - m1);
                double m3 = h * F3(x + h / 2, y1 + 2 * k2 - k1, y2 + 2 * l2 - l1, y3 + 2 * m2 - m1);

                y1 = y1 + (k1 + 4 * k2 + k3) / 6;
                y2 = y2 + (l1 + 4 * l2 + l3) / 6;
                y3 = y3 + (m1 + 4 * m2 + m3) / 6;

                x = x + h;
            }
        }

        // Функції, які описують систему диференціальних рівнянь
        static double u1(double x)
        {
            double a = 0.003;
            if (x % (2 * a) <= a)
                return 10;
            else
                return (-(x % (2 * a)) * 10 / a + 10);
        }

        static double L2(double i3)
        {
            double a0 = (L_min * i_max * i_max - L_min * i_min * i_min - L_max * i_max * i_max + L_max * i_min * i_min) / (i_max * i_max - i_min * i_min);
            double a1 = a0 + R1 * i_min;
            double a2 = a0 + R1 * i_max;
            double a3 = a0 + R2 * i_max;

            if (Math.Abs(i3) <= 1)
                return 15;
            else if (Math.Abs(i3) <= 2)
                return a0 + a1 * (Math.Abs(i3)) + a2 * (Math.Pow(Math.Abs(i3), 2)) + a3 * (Math.Pow(Math.Abs(i3), 3));
            else
                return 1.5;
        }

        static double F1(double t, double U_C1, double i3, double U_C2)
        {
            return (u1(t) - U_C1 - U_C2) * 1e3 / (R1 * C1);
        }

        static double F2(double t, double U_C1, double i3, double U_C2)
        {
            return (U_C2 - i3 * (R2 + R3)) / L2(i3);
        }

        static double F3(double t, double U_C1, double i3, double U_C2)
        {
            return (u1(t) - U_C1 - U_C2 - R1 * i3) * 1e3 / (R1 * C2);
        }

        public static void Run()
        {
            double x0 = 0;
            double[] y0 = { 1, 0, 0 };
            double xEnd = 0.03;
            double h = 0.000015;

            RungeKuttaSystem(x0, y0, xEnd, h);
        }
    }

}
