// Интерфейс квадратурных формул
trait Quadrature {
    const NAME: &'static str;

    // Порядок аппроксимации квадратурной формулы
    const ORDER: u32;

    // Вычисление интеграла
    fn integrate(&self, f: &dyn Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64;
}

// Метод Симпсона
struct Simpson;

impl Quadrature for Simpson {
    const NAME: &'static str = "Метод Симпсона";

    const ORDER: u32 = 4;

    fn integrate(&self, f: &dyn Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {

        let h = (b - a) / n as f64;

        let mut sum = h * f(a) + h * f(b);

        for i in 1..=n {
            let x_prev = a + (i - 1) as f64 * h;
            let x_curr = a + i as f64 * h;

            let midpoint = (x_prev + x_curr) / 2.0;
            sum += 4.0 * h * f(midpoint);
        }

        for i in 1..n {
            let x = a + i as f64 * h;
            sum += 2.0 * h * f(x);
        }

        sum / 6.0
    }
}

// Двухточечная квадратура Гаусса
struct Gauss2;

impl Quadrature for Gauss2 {
    const NAME: &'static str = "Гаусс-2";

    const ORDER: u32 = 4;

    fn integrate(&self, f: &dyn Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let h = (b - a) / n as f64;

    let zeta = [
        -1.0 / 3.0_f64.sqrt(),
        1.0 / 3.0_f64.sqrt(),
    ];

    let weights = [1.0, 1.0];

    let mut sum = 0.0;

    for k in 1..=n {
        let x_prev = a + (k - 1) as f64 * h;
        let h_k = h;

        let mut inner_sum = 0.0;

        for i in 0..2 {
            let x = x_prev + h_k * (zeta[i] + 1.0) / 2.0;

            inner_sum += weights[i] * f(x);
        }

        sum += h_k * inner_sum / 2.0;
    }

    sum
}
}

// Универсальный класс численного интегрирования
struct Integrator<Q: Quadrature> {
    quadrature: Q,
}

impl<Q: Quadrature> Integrator<Q> {
    fn new(quadrature: Q) -> Self {
        Self { quadrature }
    }

    fn integrate(&self, f: &dyn Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
        self.quadrature.integrate(f, a, b, n)
    }

    fn order(&self) -> u32 {
        Q::ORDER
    }

    fn name(&self) -> &'static str {
        Q::NAME
    }
}

// Аналитическое значение интеграла
fn exact_integral() -> f64 {
    0.17361304227891
}

// Результаты и таблица

fn calculate<Q: Quadrature>(integrator: &Integrator<Q>, f: &dyn Fn(f64) -> f64, a: f64, b: f64, n: usize) {
    let k = integrator.order();
    let exact = exact_integral();

    // Вложенные сетки
    let values = [
        integrator.integrate(f, a, b, n),
        integrator.integrate(f, a, b, 2 * n),
        integrator.integrate(f, a, b, 4 * n),
        integrator.integrate(f, a, b, 8 * n),
    ];

    let denominator = 2.0_f64.powi(k as i32) - 1.0;

    println!("\n{}", integrator.name());
    println!("Порядок аппроксимации k = {}", k);
    println!("Точное значение = {:.12e}\n", exact);

    println!(
        "{:>10} {:>17} {:>17} {:>25} {:>17} {:>17} {:>17}",
        "h", "I^h", "I* - I^h", "Отношение погрешностей", "Поправка", "I^R", "I* - I^R");

    // Вывод результатов
    for i in 0..3 {
        let current = values[i];

        let h = (b - a) / (n * 2_usize.pow(i as u32)) as f64;

        let error = exact - current;

        let next = values[i + 1];

        // Отношение погрешностей
        let next_error = exact - next;
        let ratio = error / next_error;

        // Поправка Ричардсона
        let correction = (next - current) / denominator;

        // Экстраполированное значение
        let richardson = next + correction;

        // Остаточная погрешность
        let richardson_error = exact - richardson;

        println!(
            "{:>10.4e} {:>17.8e} {:>17.8e} {:>25.8e} {:>17.8e} {:>17.8e} {:>17.8e}",
            h, current, error, ratio, correction, richardson, richardson_error);
    
    }
}


fn main() {
    let a = -0.07;
    let b = 0.54;
    let n = 4;

    let f = |x: f64| x.exp()-1.0;

    println!("Функция: exp(x)-1");
    println!("Интервал: [{}, {}]", a, b);
    println!("Начальное количество интервалов: {}", n);
    println!("Точное значение интеграла: {:.12e}", exact_integral());

    let simpson = Integrator::new(Simpson);

    calculate(&simpson, &f, a, b, n);

    let gauss = Integrator::new(Gauss2);

    calculate(&gauss, &f, a, b, n);
}
