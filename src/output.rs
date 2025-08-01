use core::fmt;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

use anyhow::Result;

use crate::slowdust::LCR;


pub const BUFF_SIZE: usize = 1 << 20;


pub fn write_lcr(palins: &mut Vec<LCR>, file_name: &str) -> Result<()> {
    let output = File::create(file_name)?;
    let mut writer = BufWriter::with_capacity(BUFF_SIZE, output);

    let _ = writeln!(
        writer,
        "Name\tStart\tEnd\tString\n"
    );
    for palin in palins {
        let _ = writeln!(writer, "{}", palin);
    }
    writer.flush()?;

    Ok(())
}
