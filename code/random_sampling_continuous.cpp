#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <chrono>

// Генерация простых чисел
std::vector<int> generate_primes(int limit)
{
    std::vector<bool> sieve(limit + 1, true);
    std::vector<int> primes;

    sieve[0] = sieve[1] = false;

    for (int i = 2; i <= limit; i++)
    {
        if (sieve[i])
        {
            primes.push_back(i);
            for (long long j = (long long)i * i; j <= limit && j > 0; j += i)
            {
                sieve[j] = false;
            }
        }
    }

    return primes;
}

// Структура для хранения минимальной информации о сэмпле
struct Sample
{
    int current_sum;               // Текущая сумма последовательности
    bool is_valid;                 // Валидность сэмпла
    std::vector<int> prime_values; // Значения только для простых индексов (начиная с 2)

    Sample() : current_sum(1), is_valid(true)
    {
        // x₁ = 1 по условию, но не храним его в prime_values, так как 1 не простое число
    }
};

int main()
{
    // Максимальное значение n
    const int max_n = 10000000;

    // Предварительно вычисляем простые числа
    const int prime_limit = 10000000; // Увеличиваем для надежности
    std::cout << "Generating prime numbers up to " << prime_limit << "..." << std::endl;
    auto primes = generate_primes(prime_limit);
    std::cout << "Generated " << primes.size() << " prime numbers." << std::endl;

    // Быстрая проверка на простоту
    std::vector<bool> is_prime(std::max(prime_limit, max_n) + 1, false);
    for (int p : primes)
    {
        if (p <= std::max(prime_limit, max_n))
        {
            is_prime[p] = true;
        }
    }

    // Получение индекса простого числа в массиве prime_values (начиная с 0)
    std::vector<int> prime_value_index(std::max(prime_limit, max_n) + 1, -1);
    int count = 0;
    for (int p : primes)
    {
        if (p <= std::max(prime_limit, max_n))
        {
            prime_value_index[p] = count++;
        }
    }

    // Настройка генератора случайных чисел
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);

    // Количество сэмплов
    const int num_pairs = 1000000;

    // Открываем файл для записи результатов
    std::ofstream result_file("../data/random_sampling_continuous.csv");
    if (!result_file.is_open())
    {
        std::cerr << "Error: Could not open file for writing." << std::endl;
        return 1;
    }

    // Записываем заголовки CSV
    result_file << "n,probability,1/sqrt(ln(n)),ratio" << std::endl;

    // Для отслеживания времени выполнения
    auto start_time = std::chrono::high_resolution_clock::now();

    // Создаем пары сэмплов
    std::vector<Sample> original_samples(num_pairs);

    // Обрабатываем каждое значение n последовательно
    for (int n = 1; n <= max_n; n++)
    {
        auto n_start_time = std::chrono::high_resolution_clock::now();

        // Показываем прогресс
        if (n % 50000 == 0 || n == 1 || n == 2)
        {
            std::cout << "Processing n=" << n << " (" << n << "/" << max_n << ")..." << std::endl;
        }

        // Если n - простое число (кроме 1), добавляем новое случайное значение
        if (n > 1 && is_prime[n])
        {
            int idx = prime_value_index[n];

            for (int i = 0; i < num_pairs; i++)
            {

                // Генерируем случайное значение для простого индекса
                int value = (dis(gen) == 0) ? -1 : 1;

                // Обновляем оригинальный сэмпл
                if (original_samples[i].is_valid)
                {
                    // Расширяем массив, если нужно
                    if (original_samples[i].prime_values.size() <= idx)
                    {
                        original_samples[i].prime_values.resize(idx + 1);
                    }

                    original_samples[i].prime_values[idx] = value;
                    original_samples[i].current_sum += value;
                    if (original_samples[i].current_sum < 0)
                    {
                        original_samples[i].is_valid = false;
                    }
                }
            }
        }
        // Если n - составное число (и не 1), вычисляем его значение
        else if (n > 1)
        {
            for (int i = 0; i < num_pairs; i++)
            {

                // Вычисляем значение для составного индекса
                if (original_samples[i].is_valid)
                {
                    int value = 1;
                    int num = n;

                    // Разложение на простые множители
                    for (int p : primes)
                    {
                        if (p * p > num)
                            break;

                        if (num % p == 0)
                        {
                            int power = 0;
                            while (num % p == 0)
                            {
                                num /= p;
                                power++;
                            }

                            // Если степень нечетная, умножаем на значение простого числа
                            if (power % 2 == 1)
                            {
                                int idx = prime_value_index[p];
                                value *= original_samples[i].prime_values[idx];
                            }
                        }
                    }

                    // Если осталось простое число > sqrt(n)
                    if (num > 1)
                    {
                        int idx = prime_value_index[num];
                        value *= original_samples[i].prime_values[idx];
                    }

                    original_samples[i].current_sum += value;
                    if (original_samples[i].current_sum < 0)
                    {
                        original_samples[i].is_valid = false;
                    }
                }

            }
        }

        // Подсчитываем успешные сэмплы
        int original_success = 0;

        for (int i = 0; i < num_pairs; i++)
        {
            if (original_samples[i].is_valid)
                original_success++;
        }

        // Вычисляем вероятность и теоретическое значение
        double probability = static_cast<double> (original_success) / num_pairs;
        double theoretical = (n > 1) ? 1.0 / std::sqrt(std::log(n)) : 1.0;
        double ratio = (n > 1) ? probability * std::sqrt(std::log(n)) : probability;

        // Записываем результаты в файл
        result_file << n << ","
        // << std::fixed << std::setprecision(8)
                    << probability << "," << theoretical << "," << ratio << std::endl;

        // Сбрасываем буфер файла и выводим прогресс
        if (n % 50000 == 0 || n == 1 || n == 2)
        {
            result_file.flush();

            auto n_end_time = std::chrono::high_resolution_clock::now();
            auto n_duration = std::chrono::duration_cast<std::chrono::milliseconds>(n_end_time - n_start_time).count();

            auto current_time = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();
            double avg_time_per_n = elapsed / static_cast<double>(n);
            double estimated_remaining = avg_time_per_n * (max_n - n);

            std::cout << "n=" << n << ", probability=" << probability << ", ratio=" << ratio << std::endl;
            std::cout << "Original success rate: " << static_cast<double>(original_success) / num_pairs << std::endl;
            std::cout << "Time for this n: " << n_duration / 1000.0 << " seconds" << std::endl;
            std::cout << "Elapsed time: " << elapsed / 60 << " minutes " << elapsed % 60 << " seconds" << std::endl;
            std::cout << "Estimated remaining time: " << static_cast<int>(estimated_remaining) / 60
                      << " minutes " << static_cast<int>(estimated_remaining) % 60 << " seconds" << std::endl;
            std::cout << std::endl;
        }
    }

    result_file.close();
    std::cout << "Computation complete! Results saved to random_sampling_continuous.csv" << std::endl;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    std::cout << "Total execution time: " << total_time / 60 << " minutes " << total_time % 60 << " seconds" << std::endl;

    return 0;
}
