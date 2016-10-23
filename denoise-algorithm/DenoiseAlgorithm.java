import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.badlogic.gdx.audio.analysis.FFT;

/**
 * Created by alanchern on 10/23/16.
 * Noise-reduction algorithms.
 */

public class DenoiseAlgorithm {

	public static final int SUBTRACTION = 0;
    public static final int MLSA = 1;
    public static final int MAPA = 2;
    public static final int GMAPA = 3;

    private int fftLength;
	private FFT fft;

	DenoiseAlgorithm(int fftLength, int sampleFrequency) {
		this.fftLength = fftLength;
		fft = new FFT(fftLength, sampleFrequency);
	}

	public double[] getSigma(short[] sigmaData) {
		// only need data points enough for up to five fft frames
		if (sigmaData.length > (fftLength * 3)) {
			sigmaData = Arrays.copyOfRange(sigmaData, 0, fftLength * 3);
		}
		return processSigmaData(sigmaData);
	}

	public short[] denoise(short[] data, short[] prevData, double[] sigma, int denoiseFlag) {
		// append data from previous frame, if available
		short[] inputData = append(data, prevData);

		// perform fft, process result, and peform ifft for output audio data
		short[] output = processData(inputData, prevData, sigma, denoiseFlag);
		return output;
	}

	private short[] append(short[] data, short[] prevData) {
		short[] inputData;
		if (prevData != null) {
			inputData = new short[prevData.length + data.length];
			System.arraycopy(prevData, 0, inputData, 0, prevData.length);
			System.arraycopy(data, 0, inputData, prevData.length, data.length);
			return inputData;
		} else {
			return data;
		}
	}

	private double[] processSigmaData(short[] inputData) {
		// one-time initialization of local variables
		List<float[]> fftResult;
		double[] sum = new double[fftLength];
        double[] output = null;
		int count = 0;

		// loop over data (excluding last half-frame) at half-frame intervals
		for (int i = 0; i < (inputData.length - (fftLength / 2)); i += (fftLength / 2)) {
			// check for data length to perform fft on
			if ((i + fftLength) > inputData.length) {
				fftResult = performFft(Arrays.copyOfRange(inputData, i, inputData.length));
			} else {
				fftResult = performFft(Arrays.copyOfRange(inputData, i, i + fftLength));
			}

			// get magnitude from fft result
			float[] real = fftResult.get(0);
            float[] imag = fftResult.get(1);
			double[] mag = getMagnitude(real, imag);

			for (int j = 0; j < sum.length; j++) {
				sum[j] += mag[j];
			}

			count++;
			if (count == 5) {
				output = new double[sum.length];
				for (int j = 0; j < output.length; j++) {
					output[j] = sum[j] / 5;
				}
				break;
			}
		}
		return output;
	}

	private short[] processData(short[] inputData, short[] prevData, double[] sigma, int flag) {
		// one-time initialization of local variables
		List<float[]> fftResult;
		double[] phase;
        double[] modifiedMag;
		short[] output = new short[inputData.length];
		boolean firstFrame = true;
		int count = 0;

		// loop over data (excluding last half-frame) at half-frame intervals
		for (int i = 0; i < (inputData.length - (fftLength / 2)); i += (fftLength / 2)) {
			// check for data length to perform fft on
			if ((i + fftLength) > inputData.length) {
				fftResult = performFft(Arrays.copyOfRange(inputData, i, inputData.length));
			} else {
				fftResult = performFft(Arrays.copyOfRange(inputData, i, i + fftLength));
			}

			// get phase from fft result
			float[] real = fftResult.get(0), imag = fftResult.get(1);
			double[] mag = getMagnitude(real, imag);
			phase = getPhase(real, imag);

			// de-noise algorithms based on user input
			modifiedMag = new double[mag.length];
			double gain, clean, ksi, gamma, alpha, sqrtInput;

			for (int j = 0; j < modifiedMag.length; j++) {
				if (mag[j] == 0) {
					modifiedMag[j] = 0;
				} else {
					switch (flag) {
					case SUBTRACTION:
						modifiedMag[j] = ((mag[j] - sigma[j]) < 0) ? 0 : (mag[j] - sigma[j]);
						break;
					case MLSA:
						double squaredDiff = (mag[j] * mag[j]) - (sigma[j] * sigma[j]);
						if (squaredDiff < 0) {
							squaredDiff = 0;
						}
						gain = (1 + Math.sqrt(squaredDiff / (mag[j] * mag[j]))) / 2;
						modifiedMag[j] = mag[j] * gain;
						break;
					case MAPA:
						clean = ((mag[j] - sigma[j]) < 0) ? 0 : (mag[j] - sigma[j]);
						ksi = (clean * clean) / (sigma[j] * sigma[j]);
						gamma = (mag[j] * mag[j]) / (sigma[j] * sigma[j]);
						gain = (ksi + Math.sqrt((ksi * ksi) + ((1 + ksi) * ksi / gamma))) / (2 * (1 + ksi));
						modifiedMag[j] = mag[j] * gain;
						break;
					case GMAPA:
						clean = ((mag[j] - sigma[j]) < 0) ? 0 : (mag[j] - sigma[j]);
						ksi = (clean * clean) / (sigma[j] * sigma[j]);
						gamma = (mag[j] * mag[j]) / (sigma[j] * sigma[j]);
						alpha = 2;
						sqrtInput = (ksi * ksi) + (((2 * alpha) - 1) * (alpha + ksi) * ksi / gamma);
						if (sqrtInput < 0) {
							sqrtInput = 0;
						}
						gain = (ksi + Math.sqrt(sqrtInput)) / (2 * (alpha + ksi));
						modifiedMag[j] = mag[j] * gain;
						break;
					default:
						modifiedMag[j] = 0;
						break;
					}
				}
			}

			// perform ifft on modified magnitude and unmodified phase data
			float[] ifftResult = performIfft(modifiedMag, phase);

			// check if intended output data length is greater than desired fft length
			if (output.length > fftLength) {
				// if the first frame of current iteration
				if (firstFrame) {
					firstFrame = false;
					// if no data from previous frame (very first frame)
					if (prevData == null) {
						for (int k = 0; k < ifftResult.length; k++) {
							output[k] = (short) ifftResult[k];
							count++;
						}
					} else {
						// if data from previous frame is available loop over first half of ifft result data and add on previous data
						for (int k = 0; k < (ifftResult.length / 2); k++) {
							output[k] = (short) (prevData[k] + ifftResult[k]);
							count++;
						}
						// loop over second half of ifft result data
						for (int k = (ifftResult.length / 2); k < ifftResult.length; k++) {
							output[k] = (short) ifftResult[k];
							count++;
						}
					}
				} else {
					// if not the first frame of current iteration loop over first half of ifft result data and add on previous data
					for (int k = 0; k < (ifftResult.length / 2); k++) {
						if (count < output.length) {
							output[k + (count - (ifftResult.length / 2))] += (short) ifftResult[k];
						} else {
							break;
						}
					}
					// loop over second half of ifft result data
					for (int k = (ifftResult.length / 2); k < ifftResult.length; k++) {
						if (count < output.length) {
							output[count] = (short) ifftResult[k];
							count++;
						} else {
							break;
						}
					}
				}
			} else {
				// if intended output data length is less than/equal to desired fft length, use ifft result values directly
				for (int k = 0; k < output.length; k++) {
					output[k] = (short) ifftResult[k];
				}
			}
		}
		return output;
	}

	private double[] getMagnitude(float[] real, float[] imag) {
		double[] mag = new double[real.length];
		for (int i = 0; i < mag.length; i++) {
			mag[i] = Math.sqrt((real[i] * real[i]) + (imag[i] * imag[i]));
		}
		return mag;
	}

	private double[] getPhase(float[] real, float[] imag) {
		double[] phase = new double[real.length];
		for (int i = 0; i < real.length; i++) {
			phase[i] = Math.atan2(imag[i], real[i]);
		}
		return phase;
	}

	private List<float[]> performFft(short[] data) {
		float[] modifiedData = new float[fftLength];

		// zero-pad signal
		for (int i = 0; i < fftLength; i++) {
			if (i < data.length) {
				modifiedData[i] = data[i];
			} else {
				modifiedData[i] = 0;
			}
		}

		fft.window(FFT.HAMMING);
		fft.forward(modifiedData);
		float[] real = fft.getRealPart();
		float[] imag = fft.getImaginaryPart();
		List<float[]> result = new ArrayList<>();
		result.add(real);
		result.add(imag);

        // fft result as list with first element as real part and second element as imaginary part
		return result;
	}

	private float[] performIfft(double[] mag, double[] phase) {
		float[] real = new float[mag.length], imag = new float[mag.length];
		float[] ifftResult = new float[fftLength];

		// reconstruction
		for (int i = 0; i < real.length; i++) {
			real[i] = (float) (mag[i] * Math.cos(phase[i]));
			imag[i] = (float) (mag[i] * Math.sin(phase[i]));
		}

		fft.inverse(real, imag, ifftResult);
		return ifftResult;
	}
}
