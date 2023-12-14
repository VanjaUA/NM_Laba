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
}
