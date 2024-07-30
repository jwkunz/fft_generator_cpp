#include <iostream>
#include "generated_fft.hpp"
#include <vector>
#include <chrono>
#include <random>

void naive_fft(std::vector<cfp_t> &a)
{
    const size_t N = a.size();
    if (N <= 1)
        return;

    // Divide
    std::vector<cfp_t> a_even(N / 2);
    std::vector<cfp_t> a_odd(N / 2);
    for (size_t i = 0; i < N / 2; ++i)
    {
        a_even[i] = a[i * 2];
        a_odd[i] = a[i * 2 + 1];
    }

    // Conquer
    naive_fft(a_even);
    naive_fft(a_odd);

    // Combine
    for (size_t k = 0; k < N / 2; ++k)
    {
        cfp_t t = (cfp_t)std::polar(1.0, -2 * 3.1415926535 * k / N) * a_odd[k];
        a[k] = a_even[k] + t;
        a[k + N / 2] = a_even[k] - t;
    }
}

int main(int argc, char *argv[])
{
    size_t N;

    if (argc > 1)
    {
        N = std::atoi(argv[1]);
        std::cout << "Size of the FFT to test: " << N << std::endl;
    }
    else
    {
        std::cout << "Provide the size of the FFT to test as a command line argument" << std::endl;
        return 0;
    }

    auto begin = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
    const size_t N_loops = 1000;

    // Seed for the random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-10.0, 10.0); // Range for random numbers

    std::vector<cfp_t> test_data(N);
    for (size_t i = 0; i < N; ++i)
    {
        cfp_t complexNumber(dis(gen), dis(gen));
        test_data[i] = complexNumber;
    }

    std::vector<cfp_t> naive_result(test_data);
    begin = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < N_loops; ++i)
    {
        naive_fft(naive_result);
    }
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / N_loops;
    std::cout << "Naive FFT took " << elapsed << " ns for a speed of " << 1E9 * N / elapsed << " samples per second \n";
    naive_result = test_data;
    naive_fft(naive_result);

    std::vector<cfp_t> test_result(test_data.size());
    begin = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < N_loops; ++i)
    {
        generated_fft(test_data.data(), test_result.data());
    }
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / N_loops;
    std::cout << "Generated FFT took " << elapsed << " ns for a speed of " << 1E9 * N / elapsed << " samples per second \n";

    for (size_t i = 0; i < test_data.size(); ++i)
    {
        double error = std::norm(test_result[i] - naive_result[i]);
        if (error > 1E-6)
        {
            std::cout << "Failed" << std::endl;
        }
    }
    std::cout << "Generated matched Naive FFT output" << std::endl;
}