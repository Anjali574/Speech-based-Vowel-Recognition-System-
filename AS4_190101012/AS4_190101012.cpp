// Vowel_Recognition.cpp : Defines entry point for console application.

// Included all header files required.
#include "stdafx.h"       
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>

using namespace std;

string normalized = "FALSE";                            // All global variables for storing different parameters of speech recognition.
long double samples=0, bitsPerSample=16, channel=1, sampleRate=16000;  // Speech parameters.  
int frame_size = 320;  // Frame size of frames taken.
int frame_cnt = 5; // Number of frames to be taken in one recording of vowel.
long double scaleAmplitude = 100; // Amplitude to be scaled during normalization.
long double tokhuraWt[12] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};   // Tokhura Weights.
vector<long double> raisedSineWt;    //  Raised Sine Window weights.
long double M_PI =  3.141592653589793238;   

//Get raised sine window coefficients.
void Calc_RSW(){
	raisedSineWt.clear();
	for(int i=0;i<12;i++){
		long double m = (long double)(i+1);    // Applying formula 1 + [p/2 + sin(m * PI/p)].
		long double theta = (M_PI * m)/12.0;
		long double value = sin(theta);
		value = (6.0 * value)+1.0;
		raisedSineWt.push_back(value);
	}
}

// Applying raised sine window on Ci's.
void RSW_Amplitude(vector<long double> &C){
	for(int i=1;i<=12;i++){
		C[i] = raisedSineWt[i-1]*C[i];    
	}
}

// Read contents of recordings
bool read_File(vector<long double> &amplitude, string fName){
	fstream curr;        
	curr.open(fName);    
	
    // File not found.
	if(!curr){                            
		cout<<"Failed to open file "<<fName<<"\n";
		curr.close();
		return false;
	}

	string word;
	 /*
     read file and store vowel's Samples, Amplitudes, Bits Per Sample, Channels, Sample Rate, Normalise
     */
	while(curr >> word){                   
    	if(word == "SAMPLES:"){
			curr >> word;
			samples = (long double)stoi(word);
		}
        else if(word == "CHANNELS:"){
			curr >> word;
			channel = (long double)stoi(word);
		}
		else if(word == "SAMPLERATE:"){
			curr >> word;
			sampleRate = (long double)stoi(word);
		}
		else if(word == "BITSPERSAMPLE:"){
			curr >> word;
			bitsPerSample = (long double)stoi(word);
		}
		else if(word == "NORMALIZED:"){
			curr >> word;
			normalized = word;
		}
		else{
			samples++;
			long double amp = stold(word);
			amplitude.push_back(amp);
		}
	}
	curr.close();
	return true;
}

// Write contents of reference data into file.
bool WriteFile(vector<vector<long double>> &Avg_Ci, string fName){
	ofstream Ostream(fName);
	 
	if(!Ostream){                            
		cout<<"Failed to open file "<<fName<<"\n";
		Ostream.close();
		return false;
	}

	int row = Avg_Ci.size();
	int col = Avg_Ci[0].size();
	for(int i=0;i<row;i++){
		for(int j=1;j<col;j++){
			Ostream<<Avg_Ci[i][j]<<" ";
		}
		Ostream<<"\n";
	}
	Ostream.close();
	return true;
}
//DC Shift 
void DCShift(vector<long double> &amplitude){
    long double DCOffset = 0;
	int size = amplitude.size();
	
	for(int i=0;i<size;i++)          // sum of all amplitudes 
		DCOffset += amplitude[i];
	DCOffset /= (long double)size;    //take avg of sum values.
	for(int i=0;i<size;i++){
		amplitude[i] -= DCOffset;  // Subtract DCOffset from each of amplitude values.
	}
	return;
}


// Normalise Sample Data
void Normalise(vector<long double> &amplitude){
    int size = amplitude.size();
	long double amplitudeMax = 0.0;
	
	for(int i=0;i<size;i++)
		amplitudeMax = max(amplitudeMax, amplitude[i]);   // Get max amplitude value 

	for(int i=0;i<size;i++){
		amplitude[i] /= amplitudeMax;        // Divide by max value and multiply to scale amplitude.
		amplitude[i] *= scaleAmplitude;
	}
	return;
}


// Trimming amplitude based on short term energy.
bool Trim(vector<long double> &trimAmplitude, vector<long double> &amplitude){
	long double ST_energy = 0.0;
	int size = amplitude.size(),end = frame_size*frame_cnt;
	if(size<end){        //recording very small.
		return false;
	}
	int start = 0;
	for(int i=0;i<end;i++){      // Initial recording STE.
		ST_energy += (amplitude[i]*amplitude[i]);
	}
	
	for(int i=end;i+end<size;i+=end){     // To find frame having max STE.
		long double currSTE = 0.0;
		for(int j=i;j<i+end;j++){
			currSTE += (amplitude[j]*amplitude[j]);
		}
		if(currSTE > ST_energy){
			ST_energy = currSTE;
			start = i;
		}
	}
	for(int i=0;i<end;i++){           // amplitudes of vowel region with max STE
		trimAmplitude.push_back(amplitude[start+i]);
	}
	return true;
}

//Calculating cepstral coefficients
void CalculateCIS(vector<long double> &CIs, vector<long double> &AIS){
	for(int i=1;i<=12;i++){
		long double sum = 0.0;
		for(int j=1;j<i;j++){
			sum += ((((long double)j)*CIs[j]*AIS[i-j])/((long double)i));   
		}
		long double C = AIS[i] + sum;
		CIs.push_back(C);
	}
	return;
}

//Calculating R'i(Durbin's Coef.).
void CalculateAIS(vector<long double> &AIS, vector<long double> &RIS){
	vector<long double> Energy(13, 0.0),K(13, 0.0),alpha(13, 0.0),prevAlpha(13, 0.0);
	Energy[0] = RIS[0];
	for(int i=1;i<=12;i++){
		if(i == 1){
			K[i] = RIS[i]/RIS[0];
		}
		else{
			long double sum = 0.0;
			for(int j=1;j<i;j++){           // Formula
				prevAlpha[j] = alpha[j];
				sum += (prevAlpha[j]*RIS[i-j]);
			}
			K[i] = (RIS[i] - sum)/Energy[i-1];
		}
		alpha[i] = K[i];
		for(int j=1;j<i;j++){
			alpha[j] = prevAlpha[j] - (K[i]*prevAlpha[i-j]);
		}
		Energy[i] = (1 - (K[i] * K[i])) * Energy[i-1];
	}
	for(int i=1;i<=12;i++){   // Set AIS with alpha values.
		AIS[i] = alpha[i];
	}
	return;
}

// Calc Tokhura Distance.
long double TokhuraDistance(vector<vector<long double>> &All_CIs, vector<vector<long double>> &C){
	long double tokDist = LDBL_MAX;
	int row = C.size(),col = C[0].size();
	for(int i=0;i<row;i++){
		long double distance = 0.0;
		for(int j=1;j<col;j++){
			long double diff = All_CIs[i][j] - C[i][j];    
			distance += (tokhuraWt[j-1]*(diff * diff));
		}
		if(distance < tokDist){
			tokDist = distance;
		}
	}
	return tokDist;
}

// Calc R'i values using trimmed amplitudes
vector<long double> CalculateRIS(vector<long double> &trimAmplitude, int start, int end){
	vector<long double> RIS,AIS(13, 0.0),CIs;

	long double N = (long double)frame_size;
	for(int i=0;i<=12;i++){
		long double sum = 0.0;
		for(int j=start;j+i<end;j++){
			sum += (trimAmplitude[j] * trimAmplitude[j+i]);  
		}
		sum =sum/ N;
		RIS.push_back(sum);
	}
	CalculateAIS(AIS, RIS);       // Calculate A'i.
	CIs.push_back(logl(RIS[0])); // First cepstral coefficient is logarithm of gain term of LPC model OR energy term R[0].
	CalculateCIS(CIs, AIS);
	return CIs;
}



// Function to recognize vowels into test data.
char Vowel_Recognize(vector<vector<long double>> &All_CIs, string refFolder){
	char output = '.',vowels[5] = {'a', 'e', 'i', 'o', 'u'};
	long double minTDistance = LDBL_MAX;
	vector<vector<long double>> C(5, vector<long double>(13, 0.0));
	string fName = "";
	
	for(int i=0;i<5;i++){
		fstream curr;       
		fName = refFolder + vowels[i] + ".txt";
		curr.open(fName);
		if(!curr){
			cout<<"Can't open "<<fName<<"\n";
			continue;
		}
		for(int j=0;j<5;j++){      // Get reference data for vowels[i].
			for(int k=1;k<=12;k++){
				string word;
				curr>>word;
				C[j][k] = stold(word);
			}
		}

		long double tokDist = TokhuraDistance(All_CIs, C);
		if(tokDist < minTDistance){    // Take minimum of tokhura distance of different frames.
			minTDistance = tokDist;
			output = vowels[i];
		}
		curr.close();
	}
	return output;
}

int _tmain(int argc, _TCHAR* argv[]){
    // First Get raised sine window coefficients.
	Calc_RSW(); 

	string folder = ".\\Training\\";  // Folder containing training data
	string refFolder = ".\\Vowel_Reference\\";  // Folder which will contain reference for each vowel as output.
	string testFolder = ".\\Test\\";  // Folder containing test data
	string fName = "";   // Filename of recording.
	char vowels[5] = {'a', 'e', 'i', 'o', 'u'};   
	string roll = "190101012";

    //calc for each vowel
	for(int i=0;i<5;i++){           
		vector<vector<long double>> Avg_Ci(5, vector<long double>(13, 0.0));   // Average C'i values of each frame of all recordings of each vowel
		vector<vector<long double>> All_CIs;// All C'i values of each recording and frame

		for(long long int j=1;j<=10;j++){    
			fName = folder + roll + "_" + vowels[i] + "_" + to_string(j) + ".txt";
			vector<long double> amplitude;
			if(!read_File(amplitude, fName)){
				return 0;
			}

			DCShift(amplitude);    
			Normalise(amplitude);

			vector<long double> trimAmplitude;
            // Trimming  recoding to only contain vowel and a bit of silence before & after vowel
			if(!Trim(trimAmplitude, amplitude)){    
				cout<<"Recording file very small.\n";
				return 0;
			}

			for(int k=0;k<frame_cnt;k++){         // Get C'i for every frame
				int start = (frame_size * k);
				int end = (frame_size * (k+1));
				vector<long double> CIs = CalculateRIS(trimAmplitude, start, end);
				RSW_Amplitude(CIs);     // Applying rRSW on C'i
				All_CIs.push_back(CIs);
			}
		}

		for(int j=0;j<(frame_cnt*10);j++){    // Take avg of all C'i for every samples
			for(int k=1;k<=12;k++){
				int x = j%frame_cnt;
				int y = k;
				Avg_Ci[x][y] += All_CIs[j][k]; 
			}
		}
		
		for(int j=0;j<frame_cnt;j++)
        {
			for(int k=1;k<=12;k++)
				Avg_Ci[j][k] = Avg_Ci[j][k]/10.0;
		}

		string outputfileName = refFolder + vowels[i] + ".txt";    // a.txt/e.txt .. etc
		WriteFile(Avg_Ci, outputfileName);
	}
	
	cout<<"Training done.Test Cases Output Below:\n\n\n";

	for(int i=0;i<5;i++){	          // Iterate through vowels's test data recordings.
		for(long long int j=11;j<=20;j++){
			fName = testFolder + roll + "_" + vowels[i] + "_" + to_string(j) + ".txt";
			vector<vector<long double>> All_CIs;

			vector<long double> amplitude;
            // Read recording
			if(!read_File(amplitude, fName)){      
				return 0;
			}

			DCShift(amplitude); 
			Normalise(amplitude);

			vector<long double> trimAmplitude;
            // Trimming  recoding to only contain vowel and a bit of silence before & after vowel
			if(!Trim(trimAmplitude, amplitude)){     
				cout<<"Recorded file is very small\n";
				return 0;
			}
			
			for(int k=0;k<frame_cnt;k++){      //C'i for each frame.
				int start = k*frame_size, end = (k+1)*frame_size;
				vector<long double> CIs = CalculateRIS(trimAmplitude, start, end);
                //Applying raised sine window on  C'i
				RSW_Amplitude(CIs);     
				All_CIs.push_back(CIs);
			}

            // Predicting output from reference data
			char out = Vowel_Recognize(All_CIs, refFolder);    
			if(out != '.')
			cout<<"Vowel in "<<fName<<" => "<<out<<"\n";
		}
		cout<<"\n";
	}

	return 0;
}