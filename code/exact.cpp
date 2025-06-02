#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

// Генерация простых чисел до n
std::vector<int> generate_primes(int n) {
    std::vector<bool> is_prime(n + 1, true);
    std::vector<int> primes;
    
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; i <= n; i++) {
        if (is_prime[i]) {
            primes.push_back(i);
            for (long long j = (long long)i * i; j <= n; j += i) {
                is_prime[j] = false;
            }
        }
    }
    return primes;
}

// Вычисление значения x_n для составного числа
int compute_x_value(int n, const std::vector<int>& primes, const std::vector<int>& prime_values) {
    if (n == 1) return 1; // x₁ всегда 1
    
    // Если n - простое число
    for (size_t i = 0; i < primes.size(); i++) {
        if (primes[i] == n) return prime_values[i];
    }
    
    // Если n - составное число
    int result = 1;
    int temp_n = n;
    
    for (size_t i = 0; i < primes.size(); i++) {
        int p = primes[i];
        if (p * p > temp_n) break;
        
        int count = 0;
        while (temp_n % p == 0) {
            temp_n /= p;
            count++;
        }

        if (count % 2 == 1) {
            result *= prime_values[i];
        }
    }
    
    if (temp_n > 1) {
        // Осталось простое число > sqrt(n)
        for (size_t i = 0; i < primes.size(); i++) {
            if (primes[i] == temp_n) {
                result *= prime_values[i];
                break;
            }
        }
    }
    
    return result;
}

// Рекурсивный перебор всех комбинаций
void enumerate_combinations(int index, int n, const std::vector<int>& primes, 
                           std::vector<int>& prime_values, 
                           int& valid_count, int& total_count) {
    if (index == primes.size()) {
        // Проверяем условие для всех сумм
        int sum = 0;
        bool valid = true;
        
        for (int i = 1; i <= n; i++) {
            int x_i = compute_x_value(i, primes, prime_values);
            sum += x_i;
            
            if (sum < 0) {
                valid = false;
                break;
            }
        }
        
        total_count++;
        if (valid) valid_count++;
        return;
    }
    
    // Перебираем значения для текущего простого индекса
    prime_values[index] = -1;
    enumerate_combinations(index + 1, n, primes, prime_values, valid_count, total_count);
    
    prime_values[index] = 1;
    enumerate_combinations(index + 1, n, primes, prime_values, valid_count, total_count);
}

int main() {
    const int max_n = 1000; // Ограничиваем из-за вычислительной сложности
    
    std::ofstream out("../data/exact.csv");
    out << "n,probability,1/sqrt(ln(n)),ratio" << std::endl;
    
    // Для n=1 вероятность всегда 1
    out << "1,1.00000000,inf,inf" << std::endl;
    
    for (int n = 2; n <= max_n; n++) {
        // Получаем простые числа до n
        std::vector<int> primes = generate_primes(n);
        std::vector<int> prime_values(primes.size());
        
        // Перебираем все комбинации значений для простых индексов
        int valid_count = 0, total_count = 0;
        enumerate_combinations(0, n, primes, prime_values, valid_count, total_count);
        
        double probability = (double)valid_count / total_count;
        double theoretical = 1.0 / std::sqrt(std::log(n));
        double ratio = probability * std::sqrt(std::log(n));
        
        out << n << "," << probability << "," << theoretical << "," << ratio << std::endl;
    }
    
    out.close();
    return 0;
}