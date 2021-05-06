use std::time::Instant;
use std::path::Path;
use bio::io::fasta;
use bio::data_structures::suffix_array::{lcp,suffix_array};
use std::str;
use succinct::rank::Rank9;
use bio::data_structures::rank_select::RankSelect;
use bv::BitVec;
use std::collections::{HashMap,HashSet};
mod translate;
use std::fs;

fn main() {
	let path_r = Path::new("data/reference/ref.fasta");
    let mut reader = fasta::Reader::from_file(path_r).expect("File not found");

	let mut nb_reads = 0;
	let mut nb_bases = 0;
	
	let mut cat_str = String::new();
	
	let mut ref_map = HashMap::new(); 
	
	let mut rank_bitvec: BitVec<u8> = BitVec::new();
	
	let mut counts: HashMap<String,u64> = HashMap::new();
	
	let mut ref_names: HashMap<u64,String> = HashMap::new();
	
	let start = Instant::now();
	
	//For Reference
	
	for result in reader.records() {
		let record = result.expect("Error during fasta record parsing");
		ref_map.insert(cat_str.len()+record.seq().len(),record.id().to_string());
		
		cat_str.push_str(str::from_utf8(record.seq()).expect("Error during fasta record parsing"));
		cat_str.push_str("$");
		
		for _x in 0..record.seq().len(){
			rank_bitvec.push(false);
		}
		rank_bitvec.push(true);
		counts.insert(record.id().to_string(),0);
		ref_names.insert(nb_reads,record.id().to_string());
		nb_reads += 1;
		nb_bases += record.seq().len();
	}
	
	
	let suff_arry = suffix_array(cat_str.as_bytes());
		
	
	let rs_maybe = RankSelect::new(rank_bitvec,1); 
	
	/* Print Checks
	println!("Number of reads: {}", nb_reads);
	println!("Number of bases: {}", nb_bases);
	println!("Concatenation of all Reads: \n{}", cat_str);
	println!("Len of SA: {}",suff_arry.len());
	println!("Index of SA: {}",suff_arry[126]);
	println!("Index of SA: {}",suff_arry[127]);
	println!("Rank of element 1300 {:?}",rs_maybe.rank(suff_arry[126]as u64));
	println!("Rank of element 1300 {:?}",rs_maybe.rank(suff_arry[127]as u64));
	*/
	
	//For Query
	
	let paths = fs::read_dir("data/query").unwrap();
	
	for path in paths {
		translate::p_trans(&path.unwrap().path());
    }
	
	let paths = fs::read_dir("translated").unwrap();
	
	let mut read_count = 0;
	
	for path in paths {
		//println!("Name: {}", path.unwrap().path().display());
		//let path_q = Path::new("query_test.fasta");
		reader = fasta::Reader::from_file(path.unwrap().path()).expect("File not found");
		
		let kmer = 5;
		
		for result in reader.records() {
			read_count+= 1;
			let record = result.expect("Error during fasta record parsing");
			let read = str::from_utf8(record.seq()).expect("Error during fasta record parsing");
			//println!("{}: {}",record.id(),read);
			let mut consensus: Vec<HashSet<u64>> = Vec::new();
			let mut mode_consensus: Vec<u64> = Vec::new();
			let mut cov_consensus: HashMap<u64,u64> = HashMap::new(); 
			
			let mut x = 0;
			let seqlen = record.seq().len();
			while x < (seqlen-kmer+1){
				let read_k = &read[x..(x+kmer)];
				//println!("{}",read_k);
				// Go through suffix array and find a beg and end interval
				let (beg,end) = b_search(cat_str.as_bytes(), &suff_arry, read_k.as_bytes());

				let mut mmp = 0;
				if beg <= end {
					//println!("Interval found: {}",read_k);
					//If an interval is found begin MMP
					let mut beg_p = beg;
					let mut end_p = end;
					
					//finding MMP
					//println!("SA index Match: Beg: {} & End: {}",beg,end);
					let mut done = false; 
					while !done{
						if (x+kmer+mmp+1) < read.len() {
							let temp_read = &read[x..(x+kmer+mmp+1)];
							let (temp_beg,temp_end) = b_search_mmp(cat_str.as_bytes(), &suff_arry, temp_read.as_bytes(),beg_p,end_p + 1);
							
							
							if temp_beg <= temp_end{
								beg_p = temp_beg;
								end_p = temp_end;
								mmp+=1;
							}
							else{
								done = true;
							}
						}
						else{
							done = true;
						}
					}
					//println!("MMP: {}",mmp);
					//println!("Beg prime: {} & End prime: {}", beg_p,end_p);
					let mut ref_set: HashSet<u64> = HashSet::new();
					for x in beg_p..(end_p+1){
						ref_set.insert(rs_maybe.rank(suff_arry[x] as u64).unwrap());
						mode_consensus.push(rs_maybe.rank(suff_arry[x] as u64).unwrap());
						*cov_consensus.entry(rs_maybe.rank(suff_arry[x] as u64).unwrap()).or_insert(0)  += (kmer + mmp) as u64;
					}
					consensus.push(ref_set.clone());
					//println!("{:?}",temp_cov);
				}	
				x += 1+mmp;
			}
			//Consesus machine 
			//println!("{:?}",consensus);
			if !consensus.is_empty() {
				let (intersection, others) = consensus.split_at_mut(1);
				let intersection = &mut intersection[0];
				for other in others {
					intersection.retain(|e| other.contains(e));
				}
				//println!("{:?}",intersection);
				//Using Intersection
				if !intersection.is_empty(){
					*counts.entry(ref_names.get(intersection.iter().next().unwrap()).unwrap().to_string()).or_insert(0)  += 1;
				}
				// Using Coverage
				else {
					let mut transcript = 0;
					let mut max = 0;
					for (k, v) in cov_consensus.iter() {
						if *v > max {
							transcript = *k;
							max = *v;
						}
					}
					*counts.entry(ref_names.get(&transcript).unwrap().to_string()).or_insert(0)  += 1;
				}
				/* Using Mode
				else {
					*counts.entry(ref_names.get(&mode(&mode_consensus)).unwrap().to_string()).or_insert(0)  += 1;
				}
				*/
			}
			//
			
		}
    }
	
	
	let elapsed = start.elapsed();
	println!("Millis: {} ms", elapsed.as_millis());
	println!("Reads Aligned: {}",read_count);
	println!("Counts HashMap: {:?}",counts);
}

fn mode(numbers: &[u64]) -> u64 {
    let mut occurrences = HashMap::new();

    for &value in numbers {
        *occurrences.entry(value).or_insert(0) += 1;
    }

    occurrences
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(val, _)| val)
        .expect("Cannot compute the mode of zero numbers")
}

fn b_search(s: &[u8], sa: &[usize], pat: &[u8]) -> (usize,usize) {
        
	//This section might be needed in the future if we need buckets 
	/*
	let sa = if pat.len() > 0 {
		&self.sa[self.get_bucket(pat)]
	} else {
		&self.sa[..]
	};
	*/
	b_search_mmp(s,sa,pat,0,sa.len())
}


fn b_search_mmp(s: &[u8], sa: &[usize], pat: &[u8], beg: usize, end: usize) -> (usize,usize) {
	let mut i = beg;
	let mut k = end;
	while i < k {
		let m = i + (k - i) / 2;
		if pat > &s[sa[m] as usize..] {
			i = m + 1;
		} else {
			k = m;
		}
	}

	let mut j = i;
	let mut k = end;
	while j < k {
		let m = j + (k - j) / 2;
		if s[sa[m] as usize..].starts_with(pat) {
			j = m + 1;
		} else {
			k = m;
		}
	}

	(i,j-1)
}
