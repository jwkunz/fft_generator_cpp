import numpy as np
from pathlib import Path
import os

class fft_generator_cpp:

    def configure(
            self,
            filename_out : Path,
            length_of_fft : int,
            phase_polarity : float,
            dataType = "double"
    ):
        self.filename_out = filename_out
        self.length_of_fft = length_of_fft
        self.phase_polarity = phase_polarity
        self.dataType = dataType
        self.dataTypeC = f"std::complex<{dataType}>"
        self.filestring = f"#include<complex>\n\n"
        self.filestring += f"typedef {self.dataTypeC} cfp_t;\n\n"
        self.filestring += f"void generated_fft({self.dataTypeC} in[{length_of_fft}],"
        self.filestring += f"{self.dataTypeC} out[{length_of_fft}])"+"{\n"

    def generate(self):

        input_array_list = self._make_array_list("in",self.length_of_fft)
        output_array_list = self._make_array_list("out",self.length_of_fft)
        self._recurse("in",input_array_list,output_array_list,0)

        self.filestring += "\n}"
        with (open(self.filename_out,"w") as ofs):
            ofs.write(self.filestring)
        return
    
    def _recurse(self,input_stem,input_labels,output_labels,depth):
        n = len(input_labels)
        if(n<=1): # Base case
            if(self.phase_polarity == 1):
                self.filestring += f"{output_labels[0]}={input_labels[0]}/{self.dataType}({self.length_of_fft});\n"
            else:
                self.filestring += f"{output_labels[0]}={input_labels[0]};\n"
            return output_labels
        elif(n%2 != 0): # Error case
            raise(NotImplementedError("Only supports powers of 2 length"))
        elif(n % 4 == 0) : # Radix 4
            phase_length = n//4
            results = []
            for i in range(4):
                args = self._setup_phase_args(i,input_labels,input_stem,depth,phase_length,4)
                results.append(self._recurse(args[0],args[1],args[2],args[3]))
            for j in range(phase_length):
                # Phase 0
                self.filestring += f"{output_labels[j+0*phase_length]} = "
                self.filestring += f"{results[0][j]}+"
                self.filestring += f"({self._make_polar_constant(j,n,1)})*{results[1][j]}+"
                self.filestring += f"({self._make_polar_constant(2*j,n,1)})*{results[2][j]}+"
                self.filestring += f"({self._make_polar_constant(3*j,n,1)})*{results[3][j]};\n"
                # Phase 1
                self.filestring += f"{output_labels[j+1*phase_length]} = "
                self.filestring += f"{results[0][j]}+"
                self.filestring += f"({self._make_polar_constant(j,n,-1j)})*{results[1][j]}+"
                self.filestring += f"({self._make_polar_constant(2*j,n,-1)})*{results[2][j]}+"
                self.filestring += f"({self._make_polar_constant(3*j,n,1j)})*{results[3][j]};\n"
                # Phase 2
                self.filestring += f"{output_labels[j+2*phase_length]} = "
                self.filestring += f"{results[0][j]}+"
                self.filestring += f"({self._make_polar_constant(j,n,-1)})*{results[1][j]}+"
                self.filestring += f"({self._make_polar_constant(2*j,n,1)})*{results[2][j]}+"
                self.filestring += f"({self._make_polar_constant(3*j,n,-1)})*{results[3][j]};\n"
                # Phase 3
                self.filestring += f"{output_labels[j+3*phase_length]} = "
                self.filestring += f"{results[0][j]}+"
                self.filestring += f"({self._make_polar_constant(j,n,1j)})*{results[1][j]}+"
                self.filestring += f"({self._make_polar_constant(2*j,n,-1)})*{results[2][j]}+"
                self.filestring += f"({self._make_polar_constant(3*j,n,-1j)})*{results[3][j]};\n"
            return output_labels
        else: # Radix 2
            phase_length = n//2
            # Even
            args_even = self._setup_phase_args(0,input_labels,input_stem,depth,phase_length,2)
            result_even = self._recurse(args_even[0],args_even[1],args_even[2],args_even[3])
            # Odd
            args_odd = self._setup_phase_args(1,input_labels,input_stem,depth,phase_length,2)
            result_odd = self._recurse(args_odd[0],args_odd[1],args_odd[2],args_odd[3])

            for j in range(phase_length):
                self.filestring += f"{output_labels[j]} = {result_even[j]}+({self._make_polar_constant(j,n,1)})*{result_odd[j]};\n"
                self.filestring += f"{output_labels[j+phase_length]} = {result_even[j]}-({self._make_polar_constant(j,n,1)})*{result_odd[j]};\n"
            return output_labels

    def _setup_phase_args(self,phase,input_labels,input_stem,depth,phase_length,radix):
        input_sliced = input_labels[phase::radix]
        phase_char = chr(97+phase)
        array_name_sliced = f"{input_stem}{phase_char}{depth}"
        array_string_sliced = f"{self.dataTypeC}{array_name_sliced}[{phase_length}];\n"
        self.filestring += array_string_sliced
        output_sliced = self._make_array_list(array_name_sliced,phase_length)
        next_depth = depth+1
        return (array_name_sliced,input_sliced,output_sliced,next_depth)

    def _make_array_list(self,label,count):
        result = []
        for i in range(count):
            result.append(f"{label}[{i}]")
        return result
    
    def _make_polar_constant(self,i,n,scale=1):
        value = np.exp(1j*self.phase_polarity*i*2*np.pi/n)*scale
        result = f"{self.dataTypeC}({np.real(value)},{np.imag(value)})"
        return result


if __name__ == "__main__":
    print("Unit test for fftGen")

    print("Generating the FFT")
    dut = fft_generator_cpp()
    size = 64
    dut.configure("generated_fft.hpp",size,-1,"double")
    result = dut.generate()

    print("Compiling the FFT")
    os.system("g++ -O3 fft_test_bench.cpp -o fft_test_bench.exe")

    print("Testing the FFT")
    os.system(f"./fft_test_bench.exe {size}")


    print("Passed all tests")

