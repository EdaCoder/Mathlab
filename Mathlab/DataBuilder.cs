using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Mathlab
{
    public class DataBuilder
    {
        public double[] DataLines(double fws=50) 
        {
            double fs = 5000; //sampling rate
            double fw = fws; //signal frequency
            double fn = 100; //noise frequency
            double n = 5; //number of periods to show   `
            double A = 100; //signal amplitude
            double N = 100; //noise amplitude
            int size = (int)(n * fs / fw); //sample size
            var t = Enumerable.Range(1, 12800).Select(p => p * 1 / fs).ToArray();
            var y = t.Select(p => (A * Math.Sin(2 * Math.PI * fw * p)) + (N * Math.Sin(2 * Math.PI * fn * p))).ToArray(); //Original
            return y;
        }
    }
}
