#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <unordered_map>
#include <iomanip>
#include <fstream>
#include <chrono>

// Генерируем простые числа до указанного предела
std::vector<int> generate_primes(int limit) {
    std::vector<bool> is_prime(limit + 1, true);
    std::vector<int> primes;
    
    is_prime[0] = is_prime[1] = false;
    
    for (int i = 2; i <= limit; i++) {
        if (is_prime[i]) {
            primes.push_back(i);
            for (long long j = (long long)i * i; j <= limit && j > 0; j += i) {
                is_prime[j] = false;
            }
        }
    }
    
    return primes;
}

// Функция для проверки простоты числа
bool is_prime(int n, const std::vector<int>& primes) {
    if (n <= 1) return false;
    if (n <= primes.back()) {
        return std::binary_search(primes.begin(), primes.end(), n);
    }
    
    for (int p : primes) {
        if (p * p > n) break;
        if (n % p == 0) return false;
    }
    return true;
}

// Разложение на простые множители
std::vector<std::pair<int, int>> prime_factorization(int n, const std::vector<int>& primes) {
    std::vector<std::pair<int, int>> factors;
    
    for (int p : primes) {
        if (p * p > n) break;
        
        if (n % p == 0) {
            int count = 0;
            while (n % p == 0) {
                n /= p;
                count++;
            }
            factors.push_back({p, count});
        }
    }
    
    if (n > 1) {
        factors.push_back({n, 1});
    }
    
    return factors;
}

// Проверяем последовательность
// invert_primes = true означает, что мы инвертируем все значения с простыми индексами
bool check_sequence(int n, std::mt19937& gen, std::uniform_int_distribution<>& dis, 
                    const std::vector<int>& primes, bool invert_primes = false) {
    std::unordered_map<int, int> x_values;
    
    // Функция для получения значения x_i
    std::function<int(int)> get_x = [&](int i) -> int {
        // Если значение уже вычислено, возвращаем его
        if (x_values.find(i) != x_values.end()) {
            return x_values[i];
        }

        // Если индекс равен 1, возвращаем 1
        if (i == 1) {
            int value = 1;
            x_values[i] = value;
            return value;
        }
        
        // Если индекс простой, генерируем случайное значение
        if (is_prime(i, primes)) {
            int value = (dis(gen) == 0) ? -1 : 1;
            // Если требуется инвертировать значения простых индексов, меняем знак
            if (invert_primes) {
                value = -value;
            }
            x_values[i] = value;
            return value;
        }
        
        // Для составных индексов вычисляем значение по формуле
        auto factors = prime_factorization(i, primes);
        int result = 1;
        for (const auto& factor : factors) {
            int prime = factor.first;
            int power = factor.second;
            int prime_value = get_x(prime);

            if (power > 1) {
                result = 0;
                break;
            }
            
            // Если степень = 1, то умножаем на значение простого числа
            if (power == 1) {
                result *= prime_value;
            }
        }
        
        x_values[i] = result;
        return result;
    };
    
    // Проверяем условие для каждой частичной суммы
    int sum = 0;
    for (int i = 1; i <= n; i++) {
        sum += get_x(i);
        if (sum < 0) {
            return false;
        }
    }
    
    return true;
}

// Функция для проверки пары дополняющих друг друга сэмплов
std::pair<bool, bool> check_stratified_pair(int n, std::mt19937& gen, std::uniform_int_distribution<>& dis, 
                                               const std::vector<int>& primes) {
    // Генерируем один seed для обоих вариантов, чтобы гарантировать одинаковую последовательность случайных чисел
    uint32_t seed = gen();
    std::mt19937 local_gen(seed);
    
    // Проверяем оригинальный сэмпл
    bool original_result = check_sequence(n, local_gen, dis, primes, false);
    
    // Сбрасываем генератор к тому же начальному состоянию
    local_gen.seed(seed);
    
    // Проверяем инвертированный сэмпл
    bool inverted_result = check_sequence(n, local_gen, dis, primes, true);
    
    return {original_result, inverted_result};
}

int main() {
    // Максимальное значение n
    const int max_n = 1000000;
    
    // Предварительно вычисляем простые числа
    const int prime_limit = 1000000;
    std::cout << "Generating prime numbers up to " << prime_limit << "..." << std::endl;
    auto primes = generate_primes(prime_limit);
    std::cout << "Generated " << primes.size() << " prime numbers." << std::endl;
    
    // Настройка генератора случайных чисел
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);
    
    // Количество пар сэмплов для каждого n
    const int num_pairs = 50000; // 5000 пар = 10000 сэмплов
    
    // Открываем файл для записи результатов
    std::ofstream result_file("../data/random_sampling_stratified_mu.csv");
    if (!result_file.is_open()) {
        std::cerr << "Error: Could not open file for writing." << std::endl;
        return 1;
    }
    
    // Записываем заголовки CSV
    result_file << "n,probability,1/sqrt(ln(n)),ratio,original_success,inverted_success,both_success,none_success" << std::endl;
    
    // Для отслеживания времени выполнения
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Обрабатываем каждое значение n последовательно
    for (int n = 1; n <= max_n; n++) {
        auto n_start_time = std::chrono::high_resolution_clock::now();
        
        // Показываем прогресс
        if (n % 50 == 0 || n == 1) {
            std::cout << "Processing n=" << n << " (" << n << "/" << max_n << ")..." << std::endl;
        }
        
        int original_success = 0;
        int inverted_success = 0;
        int both_success = 0;
        int none_success = 0;
        
        // Запускаем симуляции для пар сэмплов
        for (int pair = 0; pair < num_pairs; pair++) {
            auto [orig_result, inv_result] = check_stratified_pair(n, gen, dis, primes);
            
            if (orig_result) original_success++;
            if (inv_result) inverted_success++;
            
            if (orig_result && inv_result) both_success++;
            if (!orig_result && !inv_result) none_success++;
        }
        
        // Общее количество успешных сэмплов
        int total_success = original_success + inverted_success;
        
        // Вычисляем вероятность и теоретическое значение
        double probability = static_cast<double>(total_success) / (2 * num_pairs);
        double theoretical = (n > 1) ? 1.0 / std::sqrt(std::log(n)) : 1.0;
        double ratio = (n > 1) ? probability * std::sqrt(std::log(n)) : probability;
        
        // Записываем результаты в файл
        result_file << n << ","
                // << std::fixed << std::setprecision(8)
                  << probability << "," << theoretical << "," << ratio << ","
                  << static_cast<double>(original_success) / num_pairs << ","
                  << static_cast<double>(inverted_success) / num_pairs << ","
                  << static_cast<double>(both_success) / num_pairs << ","
                  << static_cast<double>(none_success) / num_pairs << std::endl;
        
        // Сбрасываем буфер файла после каждой записи
        result_file.flush();
        
        // Для каждого n показываем детальную информацию
        if (n % 50 == 0 || n == 1) {
            auto n_end_time = std::chrono::high_resolution_clock::now();
            auto n_duration = std::chrono::duration_cast<std::chrono::milliseconds>(n_end_time - n_start_time).count();
            
            auto current_time = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();
            double avg_time_per_n = elapsed / static_cast<double>(n);
            double estimated_remaining = avg_time_per_n * (max_n - n);
            
            std::cout << "  n=" << n << ", probability=" << probability << ", ratio=" << ratio << std::endl;
            std::cout << "  Original success: " << static_cast<double>(original_success) / num_pairs << std::endl;
            std::cout << "  Inverted success: " << static_cast<double>(inverted_success) / num_pairs << std::endl;
            std::cout << "  Both success: " << static_cast<double>(both_success) / num_pairs << std::endl;
            std::cout << "  None success: " << static_cast<double>(none_success) / num_pairs << std::endl;
            std::cout << "  Time for this n: " << n_duration / 1000.0 << " seconds" << std::endl;
            std::cout << "  Elapsed time: " << elapsed / 60 << " minutes " << elapsed % 60 << " seconds" << std::endl;
            std::cout << "  Estimated remaining time: " << static_cast<int>(estimated_remaining) / 60 
                      << " minutes " << static_cast<int>(estimated_remaining) % 60 << " seconds" << std::endl;
            std::cout << std::endl;
        }
    }
    
    result_file.close();
    std::cout << "Computation complete! Results saved to random_sampling_stratified_mu.csv" << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    std::cout << "Total execution time: " << total_time / 60 << " minutes " << total_time % 60 << " seconds" << std::endl;
    
    return 0;
}
