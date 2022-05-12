const L: f64 = 1.;

fn peak(x: f64) -> f64 {
    let r2 = x * x;
    let eps = 0.2;
    let eps2 = eps * eps;
    if r2 / eps2 < 1. {
        (1. - r2 / eps2).powi(4)
    } else {
        0.
    }
}

fn bosse(x: f64) -> f64 {
    peak(x - 0.5)
}

fn main() {
    let nx = 1000;

    let dx = L / nx as f64;

    let tmax = 3.;

    let cfl = 0.8;

    let dt = dx * cfl;

    let mut u_now = vec![0.; nx + 2];
    let mut unext = vec![0.; nx + 2];
    let mut uprev = vec![0.; nx + 2];

    let mut xc = vec![0.; nx + 2];

    for i in 0..nx + 2 {
        xc[i] = i as f64 * dx - dx / 2.;
    }

    for i in 0..nx + 2 {
        u_now[i] = bosse(xc[i]);
        uprev[i] = bosse(xc[i]);
    }

    let mut t = 0.;

    let b = dt / dx;
    let mut count = 0;
    let plotfreq = 100;

    while t < tmax {
        if count % plotfreq == 0 {
            plot1d(&xc, &u_now, &u_now);
        }
        count += 1;

        for i in 1..nx + 1 {
            unext[i] =
                -uprev[i] + 2. * (1. - b * b) * u_now[i] + b * b * (u_now[i - 1] + u_now[i + 1]);
        }
        //unext[0] = unext[1];
        // fixed the left 
        unext[0] = 0.;
        // free on the right
        unext[nx + 1] = unext[nx];
        //unext[nx + 1] = 0.;
        for i in 0..nx + 2 {
            uprev[i] = u_now[i];
            u_now[i] = unext[i];
        }
        t += dt;

        println!("t={}, dt={}", t, dt);
    }

    plot1d(&xc, &u_now, &u_now);
}

fn plot1d(x: &Vec<f64>, y: &Vec<f64>, z: &Vec<f64>) {
    let filename = "ploplo.dat";
    use std::fs::File;
    use std::io::BufWriter;
    use std::io::Write;
    {
        let meshfile = File::create(filename).unwrap();
        let mut meshfile = BufWriter::new(meshfile); // create a buffer for faster writes...

        x.iter()
            .zip(y.iter().zip(z.iter()))
            .for_each(|(x, (y, z))| {
                writeln!(meshfile, "{} {} {}", *x, *y, *z).unwrap();
            });
    }

    use std::process::Command;
    Command::new("python3")
        .arg("src/plot1d.py")
        .status()
        .expect("plot failed !");
}

