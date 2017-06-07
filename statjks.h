#pragma once

#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <cstring>
#include <cassert>

using namespace std;

template<class T>
vector<T> construct_jackknife_sample(vector<T>& input, int skip, int bin){
    typename vector<T>::const_iterator it = input.begin();
    vector<T> front(it, it+skip*bin);
    front.insert(front.end(), input.begin()+(skip+1)*bin, input.end());
    return front;
}

double auto_correlation(vector<double> &data);
double auto_correlation(vector<double> &data, vector<double> &series);
double mean(double* data, int num, double& average);
double mean1(double* data, int num, double& average);
double correlator(vector<double> &data, int m, double average);
array<double, 2> jackknife(vector<double> &data, int bin_size);
double jackknife_2_log_dividing(double* data1, double* data2, int num, int binSize, double& average);
double jackknife_2_power_dividing(double* data1, double* data2, int num, int binSize, double power1, double power2, double& average);
double jackknife_data_output(double* data, int num, int binSize, double& average, double* modified, int& modified_num);

double auto_correlation(vector<double> &data)
{
	double average = 0.;
	double tau_sum = 0.;

	double std = mean(data.data(), data.size(), average);

	double C0 = std * std;
	if(C0 < 1e-12){
		printf("C0 is too small:%.8e\n", C0);
		return -1.;
	}
	double Cn;
	for(long n = 0; n < data.size(); n++){
		Cn = correlator(data, n, average);
		// printf("atc: %d\tC_n: %.8f\n", n, Cn);
		if(Cn * Cn / (C0 * C0)< 1e-4) break;
		tau_sum += Cn / (2. * C0);
	}

	return tau_sum * 2.;
}

double auto_correlation(vector<double> &data, vector<double> &series)
{
	double average = 0.;
	double tauSum = 0.;

	double std = mean(data.data(), data.size(), average);

	series.clear();

	double C0 = std * std;
	double Cn;
	for(int n = 0; n < data.size() / 2; n++){
		Cn = correlator(data, n, average);
		if(Cn < 0) break;
		tauSum += Cn / (2 * C0);
		series.push_back(Cn / (2 * C0));
	}

	return tauSum * 2.;
}

double correlator(vector<double> &data, int m, double average){
	
	double sum = 0.;
	double Cm = 0.;
	
    if(m > 0){
#pragma omp parallel for reduction(+:sum)
		for(int j = 0; j < data.size() - m; j++){
            sum += (data[j + m] - average) * (data[j] - average);
        }
        Cm = (sum / (data.size() - m));
    }else{
#pragma omp parallel for reduction(+:sum)
        for(int j = - m; j < data.size(); j++){
            sum += (data[j + m] - average) * (data[j] - average);
        }
        Cm = (sum / (data.size() + m));
    }
	
	return Cm;
}

double mean(double* data, int num, double& average)
{
    double sum = 0;
    double sqrSum = 0;
    double std = 0.0;
    
    for(int i = 0; i < num; i++)
    {
        sum += data[i];
        sqrSum += data[i] * data[i];
    }
    
    average = sum / num;
    std = sqrt((sqrSum - num * average * average) / num);
    
    return std;
}

double mean1(double* data, int num, double& average)
{
    double sum = 0;
    double sqrSum = 0;
    double std = 0.0;
    
    for(int i = 0; i < num; i++)
    {
        sum += data[i];
        sqrSum += data[i] * data[i];
    }
    
    average = sum / num;
    std = sqrt((sqrSum - num * average * average) / (num - 1));
    
    return std;
}

array<double, 2> jackknife(vector<double> &data, int bin_size){
    
    double binning_sum = 0.;
   
	long binned_data_count = data.size() / bin_size;
	vector<double> binned_data(0); binned_data.reserve(binned_data_count);
	vector<double> knifed_data(0); knifed_data.reserve(binned_data_count);

	for(long i = 0; i < data.size(); i++){
		binning_sum += data[i];
		if((i + 1) % bin_size == 0){
			binned_data.push_back(binning_sum / bin_size);			
			binning_sum = 0.;
		}
	}
		
	double knifing_sum = 0.;
#pragma omp parallel for reduction(+:knifing_sum)
	for(long i = 0; i < binned_data.size(); i++){
        knifing_sum += binned_data[i];
    }
#pragma omp parallel for		
    for(long i = 0; i < binned_data.size(); i++){
        knifed_data[i] = (knifing_sum - binned_data[i]) / (binned_data_count - 1);
	}
	
    long double sum = 0.L;
    long double sigma_sum = 0.L;
	
#pragma omp parallel for reduction(+:sum)
    for(long i = 0; i < binned_data.size(); i++){
        sum += knifed_data[i];
    }
	long double average = sum / binned_data_count;
#pragma omp parallel for reduction(+:sigma_sum)
    for(long i = 0; i < binned_data.size(); i++){
        sigma_sum += (knifed_data[i] - average) * (knifed_data[i] - average);
    }   

	
	// printf("%.8f\t%.8f\n", average, sigma_sum);
	double error = sqrt(sigma_sum * (binned_data_count - 1) / binned_data_count);

	array<double, 2> ret; ret[0] = average; ret[1] = error;
	return ret;
}

double jackknife_data_output(double* data, int num, int binSize, double& average, double* modified, int& modified_num){
    double sum1 = 0.;
    
    double vbin[num / binSize];
    double vjack[num / binSize];
    
	for(int i = 0; i < num; i++){
		sum1 += data[i];
		
		if((i + 1) % binSize == 0){
			vbin[(i + 1) / binSize - 1] = sum1 / binSize;			
			sum1 = 0.;
		}
	}
		
	double sum2 = 0.;
	for (int k = 0; k < num / binSize; k++) {
        sum2 += vbin[k];
    }
		
    for (int k = 0; k < num / binSize; k++) {
        vjack[k] = (sum2 - vbin[k]) / (num / binSize - 1);
    	modified[k] = vjack[k];
	}
	
    double sum3 = 0.;
    double sqrsum = 0.;
	
    for (int k = 0; k < num / binSize; k++) {
        sum3 += vjack[k];
        sqrsum += vjack[k] * vjack[k];
    }
    
	modified_num = num / binSize;
    average = sum3 / (num / binSize);
    return sqrt((sqrsum - (num / binSize) * average * average) * ((num / binSize) - 1) / (num / binSize));
}

double jackknife_2_log_dividing(double* data1, double* data2, int num, int binSize, double& average){
    
    double sum1_1 = 0.;
	double sum1_2 = 0.;
    
    double vbin1[num / binSize];
    double vjack1[num / binSize];
	
    double vbin2[num / binSize];
    double vjack2[num / binSize];
    
	for(int i = 0; i < num; i++){
		sum1_1 += data1[i];
		
		if((i + 1) % binSize == 0){
			vbin1[(i + 1) / binSize - 1] = sum1_1 / binSize;			
			sum1_1 = 0.;
		}
	}
	
	for(int i = 0; i < num; i++){
		sum1_2 += data2[i];
		
		if((i + 1) % binSize == 0){
			vbin2[(i + 1) / binSize - 1] = sum1_2 / binSize;			
			sum1_2 = 0.;
		}
	}
		
	double sum2_1 = 0.;
	for (int k = 0; k < num / binSize; k++) {
        sum2_1 += vbin1[k];
    }
	
	double sum2_2 = 0.;
	for (int k = 0; k < num / binSize; k++) {
        sum2_2 += vbin2[k];
    }
		
    for (int k = 0; k < num / binSize; k++) {
        vjack1[k] = (sum2_1 - vbin1[k]) / (num / binSize - 1);
    }
    for (int k = 0; k < num / binSize; k++) {
        vjack2[k] = (sum2_2 - vbin2[k]) / (num / binSize - 1);
    }
	
    double sqrsum = 0.;
	
	double log_dividing_avg = log(sum2_1 / sum2_2);
	double dif;
	
    for (int k = 0; k < num / binSize; k++) {
		dif = log(vjack1[k] / vjack2[k]) - log_dividing_avg;
        sqrsum += dif * dif;
    }
    
    average = log_dividing_avg;
    return sqrt(sqrsum * (num / binSize - 1) / (num / binSize));
}

double jackknife_2_power_dividing(double* data1, double* data2, int num, int binSize, double power1, double power2, double& average){
    double sum1_1 = 0.;
	double sum1_2 = 0.;
    
    double vbin1[num / binSize];
    double vjack1[num / binSize];
	
    double vbin2[num / binSize];
    double vjack2[num / binSize];
    
	for(int i = 0; i < num; i++){
		sum1_1 += data1[i];
		
		if((i + 1) % binSize == 0){
			vbin1[(i + 1) / binSize - 1] = sum1_1 / binSize;			
			sum1_1 = 0.;
		}
	}
	
	for(int i = 0; i < num; i++){
		sum1_2 += data2[i];
		
		if((i + 1) % binSize == 0){
			vbin2[(i + 1) / binSize - 1] = sum1_2 / binSize;			
			sum1_2 = 0.;
		}
	}
		
	double sum2_1 = 0.;
	for (int k = 0; k < num / binSize; k++) {
        sum2_1 += vbin1[k];
    }
	
	double sum2_2 = 0.;
	for (int k = 0; k < num / binSize; k++) {
        sum2_2 += vbin2[k];
    }
		
    for (int k = 0; k < num / binSize; k++) {
        vjack1[k] = (sum2_1 - vbin1[k]) / (num / binSize - 1);
    }
    for (int k = 0; k < num / binSize; k++) {
        vjack2[k] = (sum2_2 - vbin2[k]) / (num / binSize - 1);
    }
	
    double sqrsum = 0.;
	
	double power_dividing_avg = pow(sum2_1 / (num / binSize), power1) / pow(sum2_2 / (num / binSize), power2);
	double dif;
	
    for (int k = 0; k < num / binSize; k++) {
		dif = pow(vjack1[k], power1) / pow(vjack2[k], power2) - power_dividing_avg;
        sqrsum += dif * dif;
    }
    
    average = power_dividing_avg;
    return sqrt(sqrsum * (num / binSize - 1) / (num / binSize));
}

string format_error(double value, double error){
	assert(error > 0.);
	char output[512];
	if(error > 10.){
		sprintf(output, "%d(%d)", int(value), int(error));
		return string(output);
	}else if(error > 1.){
		sprintf(output, "%0.1f(%0.1f)", value, error);
	}else{
		int d = int(2.-1E-12-log(error)/log(10.));
		char format[512];
		sprintf(format, "%%0.%df(%%d)", d);
		sprintf(output, format, value, int(error*pow(10.,d)));
	}
	return string(output);
}

