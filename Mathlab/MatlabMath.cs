using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Filtering;
using MathNet.Filtering.FIR;

namespace Mathlab
{
    public class MatlabMath
    {
        /// <summary>
        /// 带通
        /// </summary>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <param name="step">步进值</param>
        /// <param name="datas"></param>
        /// <returns></returns>
        public double[] Bandpass(double min, double max, int step, double[] datas)
        {
            var Bandpass = OnlineFirFilter.CreateBandpass(ImpulseResponse.Finite, 5000, min, max, step);
            return Bandpass.ProcessSamples(Bandpass.ProcessSamples(datas).Reverse().ToArray());
        }

        private double a1, a2, b0, b1, b2;
        private double x1, x2, y1, y2;
        /// <summary>
        /// 50HZ工频
        /// </summary>
        /// <param name="fs">采样率</param>
        /// <param name="bw">宽度</param>
        /// <param name="f0">频率</param>
        /// <returns></returns>
        public double[] Notch50(double[] datas, double f0 = 2, double fs = 5000, double bw = 2)
        {
            var octaves = (Math.Log(f0 / (f0 - bw / 2), 2) * 2);
            //质量因子
            var q = Math.Sqrt(Math.Pow(2, octaves)) / (Math.Pow(2, octaves) - 1);

            double omega = 2.0 * Math.PI * f0 / fs;
            double alpha = Math.Sin(omega) / (2 * q);

            b0 = 1;
            b1 = -2.0 * Math.Cos(omega);
            b2 = 1;
            a1 = b1 + alpha;
            a2 = b2 - alpha;

            List<double> res = new List<double>();
            for (int index = 0; index < datas.Length; index++)
            {
                double output = b0 * datas[index] + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
                x2 = x1;
                x1 = datas[index];
                y2 = y1;
                y1 = output;
                res.Add(output);
            }
            return res.ToArray();
        }

        /// <summary>
        /// 梳状滤波器
        /// </summary>
        /// <param name="datas">待滤波序列</param>
        /// <param name="fs">采样率</param>
        /// <param name="q">品质因素</param>
        /// <param name="f0">梳状频率</param>
        /// <returns></returns>
        public double[] Comb(double[] datas, double fs = 5000, double q = 1, double f0 = 100)
        {
            var delay = fs / f0;
            double decay = Math.Exp(-Math.PI / (q * delay));
            double feedback = (1 - decay) / (1 + decay);

            double prevSample = 0.0;

            List<double> res = new List<double>();
            for (int index = 0; index < datas.Length; index++)
            {
                var currentOutput = datas[index] + (feedback * prevSample);
                res.Add(currentOutput);

                prevSample = datas[index] + (feedback * currentOutput);

            }
            return res.ToArray();
        }

        /// <summary>
        /// 基线归零
        /// </summary>
        /// <param name="datas"></param>
        /// <returns></returns>
        public double[] BaseLineZero(double[] datas)
        {
            double temp = 0;
            foreach (double value in datas)
            {
                temp += value;
            }
            double average = temp / datas.Length;

            double[] y = new double[datas.Length];
            for (int i = 0; i < datas.Length; i++)
            {
                y[i] = datas[i] - average;
            }
            return y;
        }

        /// <summary>
        /// 潜伏期
        /// </summary>
        /// <param name="datas">数据</param>
        /// <param name="offset"></param>
        /// <param name="step"></param>
        /// <returns></returns>
        public int[] Latency(double[] datas, int step, int offset = 10)
        {
            int[] res = new int[3];
            double[] temp = new double[datas.Length];

            // Call drift function
            datas = Drift(datas, 162);

            // Find the index of maximum and minimum values in x
            int res_max = Array.IndexOf(datas, datas.Max());
            int res_min = Array.IndexOf(datas, datas.Min());

            res[0] = res_max - 1;
            res[1] = res_min - 1;

            if (res[0] > res[1])
            {
                if (res[2] > step)
                {
                    int j = 1;
                    for (int i = 0; i < datas.Length; i++)
                    {
                        temp[i] = -1 * datas[i];
                    }

                    for (int i = res[1]; i >= res[1] - step; i--)
                    {
                        j++;
                        int m = (i - 13) > 0 ? (i - 13) : 0;

                        if (temp[i] < 1 && Math.Abs(temp.Skip(m).Take(i - m + 1).Max()) < offset)
                        {
                            break;
                        }
                    }

                    res[2] = res[1] - j + 1;
                }
                else
                {
                    res[2] = 0;
                }
            }
            else
            {
                if (res[0] > step)
                {
                    int j = 1;

                    for (int i = res[0]; i >= res[0] - step; i--)
                    {
                        j++;
                        int m = (i - 13) > 0 ? (i - 13) : 0;

                        if (datas[i] < 1 && Math.Abs(datas.Skip(m).Take(i - m + 1).Max()) < offset)
                        {
                            break;
                        }
                    }

                    res[2] = res[0] - j + 1;
                }
                else
                {
                    res[2] = 0;
                }
            }

            if (res[2] < 0)
            {
                res[2] = 1;
            }

            return res;
        }
        private double[] Drift(double[] x, int i)
        {
            double[] y1 = MedianFilter(x, i);
            double[] y = new double[x.Length];

            for (int l = 0; l < x.Length; l++)
            {
                y[l] = x[l] - y1[l];
            }

            return y;
        }
        /// <summary>
        /// 中值滤波
        /// </summary>
        /// <param name="x"></param>
        /// <param name="i"></param>
        /// <returns></returns>
        private double[] MedianFilter(double[] x, int i)
        {
            double[] y = new double[x.Length];

            for (int j = 0; j < x.Length; j++)
            {
                int start = Math.Max(0, j - i);
                int end = Math.Min(x.Length - 1, j + i);
                int windowSize = end - start + 1;
                double[] window = new double[windowSize];

                Array.Copy(x, start, window, 0, windowSize);
                Array.Sort(window);

                if (windowSize % 2 == 0)
                {
                    y[j] = (window[windowSize / 2 - 1] + window[windowSize / 2]) / 2;
                }
                else
                {
                    y[j] = window[windowSize / 2];
                }
            }

            return y;
        }

        /// <summary>
        /// 美化波形
        /// </summary>
        /// <param name="datas"></param>
        /// <param name="media">中值滤波器阶数</param>
        /// <param name="highorder">高通滤波器阶数</param>
        /// <param name="highcutoff">高通截止频率</param>
        /// <param name="timedelay">滤波器初始化借用的数据</param>
        /// <returns></returns>
        public double[] Beautiful(double[] datas, int media, int highorder, double highcutoff, int timedelay)
        {
            double fs = 5000; // 采样率

            double[] input = new double[datas.Length];
            Array.Copy(datas, input, datas.Length);

            int startIndex = 7698 - timedelay;
            int endIndex = 8721 + media + timedelay;
            input = input.Skip(startIndex).Take(endIndex - startIndex).ToArray();
            for (int i = 0; i < input.Length; i++)
            {
                input[i] *= -1;
            }

            double[] y = MedianFilter(input, media);

            double[] z = y.Skip(media / 2).ToArray();

            double[] m = input.Skip(media / 2).ToArray();

            double[] n = new double[m.Length];
            for (int i = 0; i < m.Length; i++)
            {
                n[i] = m[i] - z[i];
            }

            double[] res = OnlineFirFilter.CreateLowpass(ImpulseResponse.Finite, fs, highcutoff, highorder << 1).ProcessSamples(n);
            int lStartIndex = timedelay - media / 2 + highorder / 2;
            int lEndIndex = n.Length - media - timedelay + highorder / 2 - 1;
            res = res.Skip(lStartIndex).Take(lEndIndex - lStartIndex + 1).ToArray();
            return res;
        }
    }
}
