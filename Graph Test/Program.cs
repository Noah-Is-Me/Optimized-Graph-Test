using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using System.Drawing;
using System.IO;
using System.Diagnostics;

namespace Graph_Test
{
    internal static class Program
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            //Application.Run(new Form1());
            
            Random random = new Random();
            Stopwatch stopwatch = new Stopwatch();

            int runCount = 1; // How many times to run the simulation?

            double mutationStdDev = 0.001;
            double mutationRateStdDev = 0.001; // 0.000000001

            double germlineMutationOffset = -mutationStdDev / 3;
            double somaticMutationOffset = -mutationStdDev / 3;  

            double startingGermlineMutationRate = 0.01;  // average for humans is 0.000000012
            double startingSomaticMutationRate  = 0.01;  // average for humans is 0.00000028
            // Standard is 0.0001

            int individualLength = 100; // average for humans is 3200000000
            int populationSize = 100;

            int generationMax = 1000;
            //int generationMax = 25000000;

            double chartMaxY = 0.05; //0.0000005
            double yIncIntervals = 0.05;

            bool applyDriftBarrier = false;

            double idealFitness = (double)populationSize * mutationStdDev * 1000f;
            double driftBarrierScaleFactor = 1;

            for (int runs = 0; runs < runCount; runs++)
            {
                int generationCount = 0;

                double[] fittestIndividual = new double[individualLength + 2];
                double highestFitness = 0;
                double lowestFitness = 0;

                fittestIndividual[0] = startingGermlineMutationRate;
                fittestIndividual[1] = startingSomaticMutationRate;

                bool mutationRateChanged = false;

                double[] germlineDataPoints = new double[generationMax + 1];
                double[] somaticDataPoints = new double[generationMax + 1];
                double[] fitnessDataPoints = new double[generationMax + 1];
                double[] loserDataPoints = new double[generationMax + 1];

                for (int i = 2; i < individualLength + 2; i++)
                {
                    fittestIndividual[i] = 0;
                }

                stopwatch.Start();

                while (generationCount < generationMax + 1)
                {
                    generationCount++;

                    mutationRateChanged = false;

                    germlineDataPoints[generationCount - 1] = fittestIndividual[0];
                    somaticDataPoints[generationCount - 1] = fittestIndividual[1];
                    fitnessDataPoints[generationCount - 1] = highestFitness;
                    loserDataPoints[generationCount - 1] = lowestFitness;

                    lowestFitness = highestFitness;

                    double germlineMutationRate = fittestIndividual[0];
                    double[] tempFittestIndividual = new double[individualLength + 2];
                    fittestIndividual.CopyTo(tempFittestIndividual, 0);

                    for (int i = 0; i < populationSize; i++)
                    {
                        double[] currentIndividual = new double[individualLength + 2];
                        fittestIndividual.CopyTo(currentIndividual, 0);

                        for (int j = 0; j < currentIndividual.Length; j++)
                        {
                            if (random.NextDouble() <= germlineMutationRate)
                            {
                                if (j == 0 || j == 1)
                                {

                                    double newMutationRate = normalDistribution(currentIndividual[j], mutationRateStdDev);
                                    newMutationRate = Math.Max(Math.Min(newMutationRate, 1), 0.0001);

                                    currentIndividual[j] = newMutationRate;


                                }
                                else
                                {
                                    double fitnessIncrease = normalDistribution(germlineMutationOffset, mutationStdDev);
                                    if (applyDriftBarrier) fitnessIncrease = applyDriftBarrierToFitness(fitnessIncrease, currentIndividual[j]);

                                    currentIndividual[j] += fitnessIncrease;
                                }
                            }
                        }

                        double somaticMutationRate = currentIndividual[1];

                        double[] currentSomaticIndividual = new double[individualLength + 2];

                        currentIndividual.CopyTo(currentSomaticIndividual, 0);


                        double currentFitness = 0;

                        for (int u = 0; u < currentSomaticIndividual.Length; u++)
                        {
                            if (random.NextDouble() <= somaticMutationRate)
                            {
                                if (u == 0 || u == 1)
                                {

                                    double newMutationRate = normalDistribution(currentIndividual[u], mutationRateStdDev);
                                    newMutationRate = Math.Max(Math.Min(newMutationRate, 1), 0.0001);
                                    currentSomaticIndividual[u] = newMutationRate;

                                }
                                else
                                {
                                    double fitnessIncrease = normalDistribution(somaticMutationOffset, mutationStdDev);
                                    if (applyDriftBarrier) fitnessIncrease = applyDriftBarrierToFitness(fitnessIncrease, currentIndividual[u]);

                                    currentIndividual[u] += fitnessIncrease;
                                }
                            }
                            if (u != 0 && u != 1)
                            {
                                currentFitness += currentSomaticIndividual[u];
                            }
                        }

                        if (currentFitness > highestFitness)
                        {

                            if (fittestIndividual[0] != currentIndividual[0] || fittestIndividual[1] != currentIndividual[1])
                            {
                                mutationRateChanged = true;
                            }

                            highestFitness = currentFitness;
                            tempFittestIndividual = currentIndividual;
                        }
                        else if (currentFitness < lowestFitness)
                        {
                            lowestFitness = currentFitness;
                        }
                    }


                    fittestIndividual = tempFittestIndividual;

                }

                Charting chart = new Charting();
                chart.CreateCharts(germlineDataPoints, somaticDataPoints, fitnessDataPoints, loserDataPoints, generationMax, chartMaxY, yIncIntervals, populationSize);
            }

            stopwatch.Stop();
            Console.WriteLine($"Execution Time:  {stopwatch.ElapsedMilliseconds} ms, {(int)(stopwatch.ElapsedMilliseconds / 1000)} s.");



            // Helper Functions

            double applyDriftBarrierToFitness(double fitnessIncrease, double currentFitness)
            {
                /*
                double d = currentFitness - idealFitness;
                double scale = driftBarrierScaleFactor * Math.Exp(-Math.Abs(d) / currentFitness);
                fitnessIncrease *= scale;
                */

                return fitnessIncrease;
            }

            double getFitness(double[] individual)
            {
                float fitness = 0;
                for (int i = 2; i < individual.Length; i++)
                {
                    fitness += (float)individual[i];
                }
                return (double)fitness;
            }

            double normalDistribution(double mean, double stdDev)
            {
                double u1 = 1.0 - random.NextDouble(); //uniform(0,1] random doubles
                double u2 = 1.0 - random.NextDouble();
                double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2); //random normal(0,1)
                double randNormal = mean + stdDev * randStdNormal;
                return randNormal;
            }

            double map(double n, double from1, double to1, double from2, double to2)
            {
                return (n - from1) / (to1 - from1) * (to2 - from2) + from2;
            }

        }
    }
    class Charting
    {
        public void CreateCharts(double[] germlineData, double[] somaticData, double[] fitnessData, double[] loserData, int generations, double chartMaxY, double yIncIntervals, int populationSize)
        {

            
            if (somaticData.Max() > chartMaxY || germlineData.Max() > chartMaxY) {
                chartMaxY = Math.Max(somaticData.Max(), germlineData.Max());

                double x = Math.Ceiling(chartMaxY/yIncIntervals) * yIncIntervals;

                chartMaxY = x;
            }
            

            double chartMinY = 0;

            //string genLabel = generations.ToString("N0", System.Globalization.CultureInfo.InvariantCulture);
            //string popLabel = populationSize.ToString("N0", System.Globalization.CultureInfo.InvariantCulture);

            string graphInfo = $"{generations}gen {populationSize}pop {DateTime.Now.ToString("MMM-dd-yyyy hh-mm-ss-fff tt")}";
            //string imagePath = @"C:\\Users\\Noah Sonfield\\source\\repos\\Graph Test\\Graph Test\\Graphs\\" + graphInfo + ".png";
            string baseDirectory = AppDomain.CurrentDomain.BaseDirectory;
            string projectDirectory = Directory.GetParent(Directory.GetParent(Directory.GetParent(baseDirectory).FullName).FullName).FullName;
            string graphsFolderPath = Path.Combine(projectDirectory, "Graphs");
            string imagePath = Path.Combine(graphsFolderPath, graphInfo + ".png");

            Console.WriteLine(imagePath);

            Chart Chart = new Chart();
            ChartArea CA = Chart.ChartAreas.Add("A1");

            Series S3 = Chart.Series.Add("Highest Fitness");
            Series S2 = Chart.Series.Add("Somatic Mutation Rate");
            Series S1 = Chart.Series.Add("Germline Mutation Rate");
            //Series S4 = Chart.Series.Add("Lowest Fitness");

            S1.ChartType = SeriesChartType.FastLine;
            S2.ChartType = SeriesChartType.FastLine;
            S3.ChartType = SeriesChartType.FastLine;
            //S4.ChartType = SeriesChartType.FastLine;

            S3.YAxisType = AxisType.Secondary;
            //S4.YAxisType = AxisType.Secondary;

            Chart.BackColor = Color.White;
            CA.BackColor = Color.White;
            CA.AxisX.Minimum = 0;
            CA.AxisY.Minimum = chartMinY;

            CA.AxisY2.Enabled = AxisEnabled.True;
            CA.AxisY2.Minimum = fitnessData.Min();
            CA.AxisY.Maximum = chartMaxY;
            CA.AxisX.Maximum = generations;

            int chartMaxY2 = (int)Math.Ceiling(fitnessData.Max());
            CA.AxisY2.Maximum = chartMaxY2;
            //CA.AxisY2.Interval = chartMaxY2 / 20;

            if (fitnessData.Max() == fitnessData.Min())
            {
                CA.AxisY2.Maximum = fitnessData.Max() + 1;
            }

            CA.AxisX.Interval = (double)generations/10;
            CA.AxisY.Interval = chartMaxY/20;
            CA.AxisX.Title = "Generations";
            CA.AxisY.Title = "Mutation Rates";
            CA.AxisY2.Title = "Fitness (Points)";
            CA.AxisX.TitleAlignment = StringAlignment.Center;
            CA.AxisY.TitleAlignment = StringAlignment.Center;
            CA.AxisY2.TitleAlignment = StringAlignment.Center;
            CA.AxisX.TitleFont = new Font("Ariel", 15, FontStyle.Bold);
            CA.AxisY.TitleFont = new Font("Ariel", 15, FontStyle.Bold);
            CA.AxisY2.TitleFont = new Font("Ariel", 15, FontStyle.Bold);
            Chart.Titles.Add("Fitness and the Evolution of the Germline and Somatic Mutation Rates Over Generations");
            Chart.Titles.ElementAt(0).Font = new Font("Ariel", 15, FontStyle.Bold);
            Chart.Size = new Size(1920, 1080);
            Chart.Series["Germline Mutation Rate"].BorderWidth = 2;
            Chart.Series["Somatic Mutation Rate"].BorderWidth = 2;
            Chart.Series["Highest Fitness"].BorderWidth = 4;
            //Chart.Series["Lowest Fitness"].BorderWidth = 4;
            Chart.Series["Germline Mutation Rate"].Color = Color.Blue;
            Chart.Series["Somatic Mutation Rate"].Color = Color.Red;
            Chart.Series["Highest Fitness"].Color = Color.DarkGreen;
            //Chart.Series["Lowest Fitness"].Color = Color.DarkOrange;
            Chart.AntiAliasing = AntiAliasingStyles.Graphics;
            Chart.TextAntiAliasingQuality = TextAntiAliasingQuality.High;
            for (int i=0;i<generations+1; i++)
            {
                S1.Points.AddXY(i, germlineData[i]);
                S2.Points.AddXY(i, somaticData[i]);
                S3.Points.AddXY(i, fitnessData[i]);
                //S4.Points.AddXY(i, loserData[i]);
            }

            
            Legend L = Chart.Legends.Add("L");
            L.LegendStyle = LegendStyle.Column;
            L.Title = "Legend";
            L.TitleAlignment = StringAlignment.Center;
            L.BackColor = Color.LightGray;
            L.TitleFont = new Font("Ariel", 15, FontStyle.Bold);
            L.Font = new Font("Ariel", 15, FontStyle.Bold);
            L.IsDockedInsideChartArea = true;
            L.DockedToChartArea = "A1";
            //L.Docking = Docking.Bottom;
            L.Docking = Docking.Top;
            L.Alignment = StringAlignment.Near;

            Chart.SaveImage(imagePath, ChartImageFormat.Png);
        }
    }
}
