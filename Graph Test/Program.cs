using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using System.Drawing;
using System.IO;
using System.Diagnostics;
using MathNet.Numerics.Distributions;

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

            double germlineMutationMean = -mutationStdDev / 1;
            double somaticMutationMean = -mutationStdDev / 1;

            int individualLength = 100; // average for humans is 3200000000
            int populationSize = 1000; // Or 100?

            double startingGermlineMutationRate = 0.001;  // average for humans is 0.000000012
            double startingSomaticMutationRate  = 0.001;  // average for humans is 0.00000028
            // Standard is 0.0001

            //double startingGermlineMutationRate = (0.000000012 * 3200000000) / individualLength;
            //double startingSomaticMutationRate = (0.00000028 * 3200000000) / individualLength;

            int generationMax = 10000;
            //int generationMax = 25000000;

            double chartMaxY = 0.05; //0.0000005
            double yIncIntervals = 0.05;

            bool applyDriftBarrier = true;

            double ZOfPositiveGermlineMutation = (0 - germlineMutationMean) / mutationStdDev;
            //double ZOfPositiveSomaticMutation = (0 - somaticMutationMean) / mutationStdDev;
            
            double probabilityOfPositiveGermlineMutation = (1 - Normal.CDF(germlineMutationMean, mutationStdDev,0));
            //double probabilityOfPositiveSomaticMutation = (1 - Normal.CDF(somaticMutationMean,mutationStdDev,0));

            double meanPositiveGermlineMutation = 0.003; // THIS IS COMPLETELY ARBITRARY! FIX THIS!

            //double meanPositiveGermlineMutation = germlineMutationMean + (mutationStdDev * (Normal.PDF(germlineMutationMean, mutationStdDev, ZOfPositiveGermlineMutation) / probabilityOfPositiveGermlineMutation));
            //double meanPositiveSomaticMutation = somaticMutationMean + mutationStdDev * (Normal.PDF(somaticMutationMean, mutationStdDev, ZOfPositiveSomaticMutation) / probabilityOfPositiveSomaticMutation);

            double idealFitness = (meanPositiveGermlineMutation * probabilityOfPositiveGermlineMutation / individualLength) * generationMax * Math.Log10(populationSize);

            Debug.WriteLine("germlineMutationMean: " + germlineMutationMean);
            Debug.WriteLine("ZOfPositiveGermlineMutation: " + ZOfPositiveGermlineMutation);
            Debug.WriteLine("probabilityOfPositiveGermlineMutation: " + probabilityOfPositiveGermlineMutation);
            Debug.WriteLine("meanPositiveGermlineMutation: " + meanPositiveGermlineMutation);
            Debug.WriteLine("positiveMean * positiveChance: " + meanPositiveGermlineMutation * probabilityOfPositiveGermlineMutation);
            Debug.WriteLine("Ideal Fitness: " + idealFitness);

            double driftBarrierScaleFactor = 200 / idealFitness;

            for (int runs = 0; runs < runCount; runs++)
            {
                int generationCount = 0;

                double[] fittestIndividual = new double[3];
                double highestFitness = 0;

                fittestIndividual[0] = startingGermlineMutationRate;
                fittestIndividual[1] = startingSomaticMutationRate;

                //bool mutationRateChanged = false;

                double[] germlineDataPoints = new double[generationMax + 1];
                double[] somaticDataPoints = new double[generationMax + 1];
                double[] fitnessDataPoints = new double[generationMax + 1];

                for (int i = 2; i < individualLength + 2; i++)
                {
                    fittestIndividual[i] = 0;
                }

                stopwatch.Start();

                while (generationCount < generationMax + 1)
                // Set tracepoint here to track the generations while running.
                // Conditional: generationCount%10000==0    ;  Message:  Generations: {generationCount}
                // NOTE: Tracepoints drastically increase runtime!
                {

                    if (generationCount % 1000 == 0) { Debug.WriteLine("Generation: " + generationCount); }

                    generationCount++;

                    

                    //mutationRateChanged = false;

                    germlineDataPoints[generationCount - 1] = fittestIndividual[0];
                    somaticDataPoints[generationCount - 1] = fittestIndividual[1];
                    fitnessDataPoints[generationCount - 1] = highestFitness;

                    highestFitness = double.NegativeInfinity;

                    double germlineMutationRate = fittestIndividual[0];
                    double[] tempFittestIndividual = new double[3];
                    fittestIndividual.CopyTo(tempFittestIndividual, 0);

                    for (int i = 0; i < populationSize; i++)
                    {
                        double[] currentIndividual = new double[3];
                        fittestIndividual.CopyTo(currentIndividual, 0);

                        for (int j = 0; j < currentIndividual.Length; j++)
                        {
                            if (j == 0 || j == 1)
                            {
                                if (random.NextDouble() <= germlineMutationRate)
                                {
                                    double newMutationRate = normalDistribution(currentIndividual[j], mutationRateStdDev);
                                    newMutationRate = Math.Max(Math.Min(newMutationRate, 1), 0.0001);

                                    currentIndividual[j] = newMutationRate;
                                }
                            }
                            else
                            {
                                double fitnessIncrease = 0;
                                for (int o=0; o<individualLength; o++) {
                                    fitnessIncrease += normalDistribution(germlineMutationMean, mutationStdDev);
                                }
                                if (applyDriftBarrier) fitnessIncrease = applyDriftBarrierToFitness(fitnessIncrease, currentIndividual[j]);
                                currentIndividual[j] += fitnessIncrease;
                            }
                        }

                        double somaticMutationRate = currentIndividual[1];

                        double[] currentSomaticIndividual = new double[3];

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
                                    double fitnessIncrease = normalDistribution(somaticMutationMean, mutationStdDev);
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

                            /*
                            if (fittestIndividual[0] != currentIndividual[0] || fittestIndividual[1] != currentIndividual[1])
                            {
                                mutationRateChanged = true;
                            }
                            */

                            highestFitness = currentFitness;
                            tempFittestIndividual = currentIndividual;
                        }
                    }


                    fittestIndividual = tempFittestIndividual;

                }

                Charting chart = new Charting();
                chart.CreateCharts(germlineDataPoints, somaticDataPoints, fitnessDataPoints, generationMax, chartMaxY, yIncIntervals, populationSize, applyDriftBarrier, idealFitness);
            }

            stopwatch.Stop();
            Debug.WriteLine($"Execution Time:  {stopwatch.ElapsedMilliseconds} ms, {(int)(stopwatch.ElapsedMilliseconds / 1000)} s.");



            // Helper Functions

            double applyDriftBarrierToFitness(double fitnessIncrease, double currentFitness)
            {
                double d = currentFitness - idealFitness;

                double tanhFunction = (-50 * Math.Tanh( (driftBarrierScaleFactor*d) + 2 )) + 50;

                //Debug.WriteLine("Scale: " + tanhFunction + " ;   Distance: " + d);

                fitnessIncrease = fitnessIncrease * tanhFunction / 100;

                return fitnessIncrease;
            }

            double normalDistribution(double mean, double stdDev)
            {
                double u1 = 1.0 - random.NextDouble(); //uniform(0,1] random doubles
                double u2 = 1.0 - random.NextDouble();
                double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2); //random normal(0,1)
                double randNormal = mean + stdDev * randStdNormal;
                return randNormal;
            }

            /*
            double getFitness(double[] individual)
            {
                float fitness = 0;
                for (int i = 2; i < individual.Length; i++)
                {
                    fitness += (float)individual[i];
                }
                return (double)fitness;
            }
            */

            /*
            double map(double n, double from1, double to1, double from2, double to2)
            {
                return (n - from1) / (to1 - from1) * (to2 - from2) + from2;
            }
            */

        }
    }
    class Charting
    {
        public void CreateCharts(double[] germlineData, double[] somaticData, double[] fitnessData, int generations, double chartMaxY, double yIncIntervals, int populationSize, bool applyDriftBarrier, double idealFitness)
        {

            int fitnessMultiplier = 1;
            
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

            Debug.WriteLine("Image Path: " + imagePath);

            Chart Chart = new Chart();
            ChartArea CA = Chart.ChartAreas.Add("A1");

            Series fitnessSeries = Chart.Series.Add("Highest Fitness");
            Series somaticSeries = Chart.Series.Add("Somatic Mutation Rate");
            Series germlineSeries = Chart.Series.Add("Germline Mutation Rate");
            //Series driftBarrierSeries = Chart.Series.Add("Ideal Fitness");

            fitnessSeries.ChartType = SeriesChartType.FastLine;
            somaticSeries.ChartType = SeriesChartType.FastLine;
            germlineSeries.ChartType = SeriesChartType.FastLine;
            //driftBarrierSeries.ChartType = SeriesChartType.FastLine;

            fitnessSeries.YAxisType = AxisType.Secondary;
            //S4.YAxisType = AxisType.Secondary;

            Chart.BackColor = Color.White;
            CA.BackColor = Color.White;
            CA.AxisX.Minimum = 0;
            CA.AxisY.Minimum = chartMinY;

            CA.AxisY2.Enabled = AxisEnabled.True;
            CA.AxisY2.Minimum = fitnessData.Min();
            CA.AxisY.Maximum = chartMaxY;
            CA.AxisX.Maximum = generations;

            int chartMaxY2 = (int)Math.Ceiling(fitnessData.Max()*fitnessMultiplier);
            CA.AxisY2.Maximum = chartMaxY2;
            //CA.AxisY2.Interval = chartMaxY2 / 20;

            if (fitnessData.Max() == fitnessData.Min())
            {
                CA.AxisY2.Maximum = fitnessData.Max()* fitnessMultiplier + 1*fitnessMultiplier;
            }

            
            // if (applyDriftBarrier) CA.AxisY2.Maximum = idealFitness + (idealFitness / 100);

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
            //Chart.Series["Ideal Fitness"].BorderWidth = 4;

            Chart.Series["Germline Mutation Rate"].Color = Color.Blue;
            Chart.Series["Somatic Mutation Rate"].Color = Color.Red;
            Chart.Series["Highest Fitness"].Color = Color.DarkGreen;
            //Chart.Series["Ideal Fitness"].Color = Color.DarkOrange;

            Chart.AntiAliasing = AntiAliasingStyles.Graphics;
            Chart.TextAntiAliasingQuality = TextAntiAliasingQuality.High;

            for (int i=0;i<generations+1; i++)
            {
                germlineSeries.Points.AddXY(i, germlineData[i]);
                somaticSeries.Points.AddXY(i, somaticData[i]);
                fitnessSeries.Points.AddXY(i, fitnessData[i]* fitnessMultiplier);

                // if (applyDriftBarrier) driftBarrierSeries.Points.AddXY(i, idealFitness);
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
