const C: f64 = -1.;

const L: f64 = 1.;

fn peak(x: f64) -> f64 {
    let r2 = x * x;
    let eps = 0.1;
    let eps2 = eps * eps;
    if r2 / eps2 < 1. {
        (1. - r2 / eps2).powi(4)
    } else {
        0.
    }
}

fn exact_sol(x: f64, t: f64) -> f64 {
    peak(x - C * t - 0.8)
}

fn plot1d(x: &Vec<f64>, y: &Vec<f64>, z: &Vec<f64>) {
    let filename = "ploplo.dat";
        use std::fs::File;
        use std::io::BufWriter;
        use std::io::{Write};    
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


fn main()  {
    let nx = 1000;

    let dx = L / nx as f64;

    let mut un = vec![0.; nx + 1];

    let mut xc = vec![0.; nx + 1];

    for i in 0..nx + 1 {
        xc[i] = i as f64 * dx;
    }

    for i in 0..nx + 1 {
        un[i] = exact_sol(xc[i], 0.);
    }

    let mut t = 0.;

    plot1d(&xc, &un, &un);

    let tmax = 0.3;

    let cfl = 0.8;

    let dt = dx / C.abs() * cfl;


    while t < tmax {
        for i in 0..nx {
            un[i] = un[i] - C * dt / dx * (un[i + 1] - un[i]);
        }
        t +=  dt;
        un[nx] = exact_sol(xc[nx], t);

        println!("t={}, dt={}", t, dt);
    }


    plot1d(&xc, &un, &un);

}
