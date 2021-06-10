use std::time::Instant;
use std::path::Path;
use bio::io::{fasta,fastq};
use bio::data_structures::suffix_array::{lcp,suffix_array};
use std::str;
use bio::data_structures::rank_select::RankSelect;
use bv::BitVec;
use std::collections::{HashMap,HashSet};
use std::fs;
use bio::alphabets::rna::revcomp;
use protein_translate::translate;
use std::collections::LinkedList;
use std::process::{Command,Stdio};

fn main() {
	let path_r = Path::new("data/ref.fasta");
    let mut reader = fasta::Reader::from_file(path_r).expect("File not found");

	let mut nb_reads = 0;
	let mut nb_bases = 0;
	
	let mut cat_str = String::new();
	
	let mut ref_map = HashMap::new(); 
	
	let mut rank_bitvec: BitVec<u8> = BitVec::new();
	
	let mut counts: HashMap<String,u64> = HashMap::new();
	
	let mut ref_names: HashMap<u64,String> = HashMap::new();
	let mut read_alignments: HashMap<String,u64> = HashMap::new();
	


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
	let path_q = Path::new("data/sample_01_100000_interleaved.fastq");
    let reader_q =  fastq::Reader::from_file(path_q).expect("File not found");
	//let paths = fs::read_dir(r"D:/data/reads_20mil_1txpg_protCoding.rep3/query").unwrap();
	let mut records = reader_q.records();
	
    while let (Some(Ok(pair1)),Some(Ok(pair2))) = (records.next(),records.next()){
		//println!("HERE IS READ: {} {}", pair1.id(),pair2.id());
		let temp_cat_str = cat_str.as_bytes();
		let p1_frame_con = vec![palign(temp_cat_str, &suff_arry, translate(pair1.seq()).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&pair1.seq()[1..]).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&pair1.seq()[2..]).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&revcomp(pair1.seq())).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&revcomp(pair1.seq())[1..]).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&revcomp(pair1.seq())[2..]).as_bytes(), &rs_maybe)];
		//println!("{:?}",p1_frame_con);

		let p2_frame_con = vec![palign(temp_cat_str, &suff_arry, translate(pair2.seq()).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&pair2.seq()[1..]).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&pair2.seq()[2..]).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&revcomp(pair2.seq())).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&revcomp(pair2.seq())[1..]).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&revcomp(pair2.seq())[2..]).as_bytes(), &rs_maybe)];
		//println!("{:?}",p2_frame_con);

		let p1_frame = bestInFrame(&p1_frame_con);
		let p2_frame = bestInFrame(&p2_frame_con);

		//println!("{} Aligns to {:?}",pair1.id(), best_in_pair(&p1_frame_con, &p2_frame_con, p1_frame, p2_frame));

	} 

	//println!("{:?}", read_alignments);
	let elapsed = start.elapsed();
	println!("Millis: {} ms", elapsed.as_millis());
	println!("Reads Aligned: {}",counts.len());
}

fn best_in_pair(p1: &Vec<Vec<(u64,f64)>>, p2: &Vec<Vec<(u64,f64)>>, frame1: usize, frame2: usize) -> Vec<u64>{
	let mut matches: Vec<u64> = Vec::new();

	if frame1 < 10 && frame2 < 10{
		for x in &p1[frame1]{
			for y in &p2[frame2]{
				if x.0 == y.0{
					matches.push(x.0);
				}
			}
		}
	}
	else if frame2 < 10 {
		for x in &p2[frame2]{
			matches.push(x.0);
		}
	}
	else if frame1 < 10 {
		for x in &p1[frame1]{
			matches.push(x.0);
		}
	}

}

fn best_in_frame(frames: &Vec<Vec<(u64,f64)>>) -> usize{
	let mut best = 10;
	let mut max = 0;

	for (i,info) in frames.iter().enumerate(){
		if info.len() > max{
			max = info.len();
			best = i;
		}
	}

	best
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

fn palign(cat_str: &[u8], suff_arry: &[usize], read:&[u8], rs_maybe: &bio::data_structures::rank_select::RankSelect) -> Vec<(u64,f64)>{
	let mut trans: Vec<(u64,f64)> = Vec::new();
	let kmer = 5;
	let mut cov_consensus: HashMap<u64,f64> = HashMap::new(); 
	let threshold:f64 = 33.3;
	
	let mut x = 0;
	let seqlen = read.len();
	while x < (seqlen-kmer+1){
		let read_k = &read[x..(x+kmer)];
		//println!("{}",read_k);
		// Go through suffix array and find a beg and end interval
		let (beg,end) = b_search(cat_str, &suff_arry, read_k);

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
					let (temp_beg,temp_end) = b_search_mmp(cat_str, &suff_arry, temp_read,beg_p,end_p + 1);
					
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
			//println!("{} {} {}",seqlen, kmer, mmp);
			for x in beg_p..(end_p+1){
				//consensus.push(HashSet::new().insert(rs_maybe.rank(suff_arry[x] as u64).unwrap()));
				//mode_consensus.push(rs_maybe.rank(suff_arry[x] as u64).unwrap());
				*cov_consensus.entry(rs_maybe.rank(suff_arry[x] as u64).unwrap()).or_insert((kmer + mmp) as f64 / seqlen as f64 * 100.0) += 0.0 ;
			}
			
		}	
		x += 1+mmp;
	}
	//Consesus machine 
	/*
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
	*/
		// Using Coverage
		//else {
	if !cov_consensus.is_empty(){
		for (protein, coverage) in cov_consensus.iter() {
			if *coverage > threshold {
				palign_result.push((*protein,*coverage));
			}
		}
		//
	}
		//}
		/* Using Mode
		else {
			*counts.entry(ref_names.get(&mode(&mode_consensus)).unwrap().to_string()).or_insert(0)  += 1;
		}
		*/
	palign_result
}
