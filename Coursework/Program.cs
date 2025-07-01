using CourseworkLibrary;

public class Coursework
{
    public static void Main(string[] args)
    {
        int nx, ny, nz;
        double lx, ly, lz, tMax, dt, alpha;
        Func<double, double, double, double> initialCondition;

        Console.WriteLine("1 - вручную\n2 - автоматически");
        int swCheck = Convert.ToInt16(Console.ReadLine());
        switch (swCheck)
        {
            case 1:
                StartWrite.Write();

                Console.WriteLine("Введите nx:");
                nx = Convert.ToInt32(Console.ReadLine());
                StartWrite.Write(nx);

                Console.WriteLine("Введите ny:");
                ny = Convert.ToInt32(Console.ReadLine());
                StartWrite.Write(nx, ny);

                Console.WriteLine("Введите nz:");
                nz = Convert.ToInt32(Console.ReadLine());
                StartWrite.Write(nx, ny, nz);

                Console.WriteLine("Введите lx:");
                lx = Convert.ToDouble(Console.ReadLine());
                StartWrite.Write(nx, ny, nz, lx);

                Console.WriteLine("Введите ly:");
                ly = Convert.ToDouble(Console.ReadLine());
                StartWrite.Write(nx, ny, nz, lx, ly);

                Console.WriteLine("Введите lz:");
                lz = Convert.ToDouble(Console.ReadLine());
                StartWrite.Write(nx, ny, nz, lx, ly, lz);

                Console.WriteLine("Введите tMax:");
                tMax = Convert.ToDouble(Console.ReadLine());
                StartWrite.Write(nx, ny, nz, lx, ly, lz, tMax);

                Console.WriteLine("Введите dt:");
                dt = Convert.ToDouble(Console.ReadLine());
                StartWrite.Write(nx, ny, nz, lx, ly, lz, tMax, dt);

                Console.WriteLine("Введите alpha:");
                alpha = Convert.ToDouble(Console.ReadLine());
                StartWrite.Write(nx, ny, nz, lx, ly, lz, tMax, dt, alpha);

                var (function, description) = StartWrite.InitialCondition(lx, ly, lz);
                initialCondition = function;
                StartWrite.Write(nx, ny, nz, lx, ly, lz, tMax, dt, alpha, description);


                break;

            default:
                nx = 10;
                ny = 10;
                nz = 10;
                lx = 1.0;
                ly = 1.0;
                lz = 1.0;
                tMax = 0.1;
                dt = 0.001;
                alpha = 0.1;
                initialCondition = (x, y, z) => 0.0;

                StartWrite.WriteDefault(nx, ny, nz, lx, ly, lz, tMax, dt, alpha);
                break;
        }

        // Граничные условия (Дирихле, Неймана и Робин)
        var boundaryConditions = new Dictionary<string, (int type, double value)>
        {
            { "x0", (0, 15.0) },       // Дирихле: u(0, y, z, t) = 20
            { "x1", (1, 0.0) },        // Неймана: du/dx(lx, y, z, t) = 0
            { "y0", (0, 30.0) },       // Дирихле: u(x, 0, z, t) = 30
            { "y1", (1, 0.0) },        // Неймана: du/dy(x, ly, z, t) = 0
            { "z0", (0, 60.0) },       // Дирихле: u(x, y, 0, t) = 40
            { "z1", (1, 0.0) }         // Неймана: du/dz(x, y, lz, t) = 0
        };

        // Коэффициенты для граничного условия Робина
        var robinCoefficients = new Dictionary<string, (double alpha, double beta)>
        {
            { "z0", (1.0, 1.0) }
        };

        List<(double x, double y, double z, double Q)> pointSources = new List<(double x, double y, double z, double Q)>
        {
            (lx / 4, ly / 4, lz / 4, 50),
            (3 * lx / 4, 3 * ly / 4, 3 * lz / 4, 50)
        };

        HeatSourceFunction heatSource = (x, y, z, t) =>
        {
            return 100 * (x * (lx - x) * y * (ly - y) * z * (lz - z)) * Math.Exp(-t);
        };

        double[,,] solution = Solve.SolveGeneral(nx, ny, nz, lx, ly, lz, tMax, dt, alpha, initialCondition, boundaryConditions, robinCoefficients, heatSource, pointSources);


        Console.WriteLine("\nВыберите способ вывода результатов:\n1) Вывод в .txt файл\n2) Визуализировать в консоли");
        string writeCheck = Console.ReadLine();
        switch (writeCheck)
        {
            default:
                Results.WriteToFile(solution);
                break;

            case "2":
                while (true)
                {
                    Console.Clear();
                    Console.WriteLine("\n1) Выбрать срез\n2) Выйти");
                    int wCheck = Convert.ToInt16(Console.ReadLine());
                    if (wCheck == 1)
                    {
                        Results.PrintSlise(solution);
                    }
                    else
                    {
                        break;
                    }
                }
                break;
        }
    }
}
