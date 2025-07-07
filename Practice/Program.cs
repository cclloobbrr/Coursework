using PracticeLibrary;

public class Practice
{
    public static void Main(string[] args)
    {
        Console.WriteLine("Добро пожаловать");

        // Параметры сетки
        int nx = Config.nx, ny = Config.ny, nz = Config.nz;
        double lx = 1.0, ly = 1.0, lz = 1.0;

        // Параметры времени
        double tMax = 0.1, dt = 0.001;

        // Коэффициенты диффузии
        double algaeDiffusion = 0.1;  // для водорослей
        double shrimpDiffusion = 0.1; // для креветок

        Func<double, double, double, double> initialAlgae = (x, y, z) => x >= 0.5 ? 1.0 : 0.0;
        //Func<double, double, double, double> initialAlgae = (x, y, z) => 1;

        Func<double, double, double, double> initialShrimps = (x, y, z) => x < 0.5 ? 1.0 : 0.0;


        // Граничные условия для водорослей (Algae)
        var algaeBoundaryConditions = new Dictionary<string, (int type, double value)>
        {
            { "x0", (1, 0.0) },
            { "x1", (1, 0.0) },
            { "y0", (1, 0.0) },
            { "y1", (1, 0.0) },
            { "z0", (1, 0.0) },
            { "z1", (1, 0.0) }
        };

        // Граничные условия для креветок (Shrimps)
        var shrimpsBoundaryConditions = new Dictionary<string, (int type, double value)>
        {
            { "x0", (1, 0.0) },
            { "x1", (1, 0.0) },
            { "y0", (1, 0.0) },
            { "y1", (1, 0.0) },
            { "z0", (1, 0.0) },
            { "z1", (1, 0.0) }
        };

        // Коэффициенты для граничного условия Робина
        var robinCoefficients = new Dictionary<string, (double alpha, double beta)>
        {
        };

        double[,,] algaeSolution;
        double[,,] shrimpsSolution;
        //algaeSolution = Solve.SolveGeneral(nx, ny, nz, lx, ly, lz, tMax, dt, algaeDiffusion, initialAlgae, algaeBoundaryConditions, robinCoefficients);
        //shrimpsSolution = Solve.SolveGeneral(nx, ny, nz, lx, ly, lz, tMax, dt, shrimpDiffusion, initialShrimps, shrimpsBoundaryConditions, robinCoefficients);

        (algaeSolution, shrimpsSolution) = Solve_2.SolveGeneral(
            nx,
            ny,
            nz,
            lx,
            ly,
            lz,
            tMax,
            dt,
            algaeDiffusion,
            shrimpDiffusion,
            initialAlgae,
            initialShrimps,
            algaeBoundaryConditions,
            shrimpsBoundaryConditions,
            robinCoefficients);


        Console.WriteLine("\nВыберите способ вывода результатов:\n1) Вывод в .txt файл\n2) Визуализировать в консоли водоросли\n2) Визуализировать в консоли креветки");
        string writeCheck = Console.ReadLine();
        switch (writeCheck)
        {
            default:
                //Results.WriteToFile(solution);
                break;

            case "2":
                while (true)
                {
                    Console.Clear();
                    Console.WriteLine("\n1) Выбрать срез\n2) Выйти");
                    int wCheck = Convert.ToInt16(Console.ReadLine());
                    if (wCheck == 1)
                    {
                        Results.PrintSlice(algaeSolution);
                    }
                    else
                    {
                        break;
                    }
                }
                break;

            case "3":
                while (true)
                {
                    Console.Clear();
                    Console.WriteLine("\n1) Выбрать срез\n2) Выйти");
                    int wCheck = Convert.ToInt16(Console.ReadLine());
                    if (wCheck == 1)
                    {
                        //Results.PrintSlice(shrimpsSolution);
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
