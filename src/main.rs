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
use argparse::{ArgumentParser, StoreTrue, Store};

fn main() {

	let mut ref_file = String::from("data/ref.fasta");
	let mut query_file = String::from("data/sample_01_100000_interleaved.fastq");
	{ // this block limits scope of borrows by ap.refer() method
	let mut ap = ArgumentParser::new();
	ap.set_description("Greet somebody.");
	ap.refer(&mut ref_file)
	.add_argument("Reference Sequence", Store, "Reference Fasta file");
	ap.refer(&mut query_file)
	.add_argument("Query Sequence", Store, "Query Interleaved Fastq file");
	ap.parse_args_or_exit();
	}
	let path_r = Path::new(&ref_file);
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
	let start_ref = Instant::now();
	let start_ref_cat = Instant::now();
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
	let elapsed_ref_cat = start_ref_cat.elapsed();
	
	let start_ref_suff = Instant::now();
	let suff_arry = suffix_array(cat_str.as_bytes());
	let elapsed_ref_suff = start_ref_suff.elapsed();
	
	let start_ref_rank = Instant::now();
	let rs_maybe = RankSelect::new(rank_bitvec,1); 
	let elapsed_ref_rank = start_ref_rank.elapsed();

	let elapsed_ref = start_ref.elapsed();
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
	let start_query = Instant::now();
	let path_q = Path::new(&query_file);
    let reader_q =  fastq::Reader::from_file(path_q).expect("File not found");
	//let paths = fs::read_dir(r"D:/data/reads_20mil_1txpg_protCoding.rep3/query").unwrap();
	let mut records = reader_q.records();
	
    while let (Some(Ok(pair1)),Some(Ok(pair2))) = (records.next(),records.next()){
		//println!("HERE IS READ: {} {}", pair1.id(),pair2.id());
		let temp_cat_str = cat_str.as_bytes();
		let start_query_frame1 = Instant::now();
		let p1_frame_con = vec![palign(temp_cat_str, &suff_arry, translate(pair1.seq()).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&pair1.seq()[1..]).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&pair1.seq()[2..]).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&revcomp(pair1.seq())).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&revcomp(pair1.seq())[1..]).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&revcomp(pair1.seq())[2..]).as_bytes(), &rs_maybe)];
		//println!("{:?}",p1_frame_con);
		let elapsed_query_frame1 = start_query_frame1.elapsed();
		//println!("Query Frame1 time: {} ms", elapsed_query_frame1.as_millis());

		let start_query_frame2 = Instant::now();
		let p2_frame_con = vec![palign(temp_cat_str, &suff_arry, translate(pair2.seq()).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&pair2.seq()[1..]).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&pair2.seq()[2..]).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&revcomp(pair2.seq())).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&revcomp(pair2.seq())[1..]).as_bytes(), &rs_maybe),
								palign(temp_cat_str, &suff_arry, translate(&revcomp(pair2.seq())[2..]).as_bytes(), &rs_maybe)];
		//println!("{:?}",p2_frame_con);
		let elapsed_query_frame2 = start_query_frame2.elapsed();
		//println!("Query Frame2 time: {} ms", elapsed_query_frame2.as_millis());

		let start_query_best_frame = Instant::now();
		let p1_frame = best_in_frame(&p1_frame_con);
		let p2_frame = best_in_frame(&p2_frame_con);

		let elapsed_query_best_frame = start_query_best_frame.elapsed();
		//println!("Query Best Frame time: {} ms", elapsed_query_best_frame.as_millis());

		let start_query_best = Instant::now();
		best_in_pair(&p1_frame_con, &p2_frame_con, p1_frame, p2_frame);
		//println!("{} Aligns to {:?}",pair1.id(), );
		let elapsed_query_best = start_query_best.elapsed();
		//println!("Query Best time: {} ms", elapsed_query_best.as_millis());
	} 
	let elapsed_query = start_query.elapsed();
	//println!("{:?}", read_alignments);
	let elapsed = start.elapsed();
	println!("time: {} seconds", elapsed.as_millis()/1000);
	println!("Read and Cat time: {} seconds", elapsed_ref_cat.as_millis()/1000);
	println!("Suff time: {} seconds", elapsed_ref_suff.as_millis()/1000);
	println!("Rank time: {} seconds", elapsed_ref_rank.as_millis()/1000);
	println!("Query time: {} seconds", elapsed_query.as_millis()/1000);
	println!("Entire Ref time: {} seconds", elapsed_ref.as_millis()/1000);
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

	matches
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
	//let start_align = Instant::now();
	let mut x = 0;
	let seqlen = read.len();
	while x < (seqlen-kmer+1){
		let read_k = &read[x..(x+kmer)];
		//println!("{}",read_k);
		// Go through suffix array and find a beg and end interval
		//let start_align_initial = Instant::now();
		let (beg,end) = b_search(cat_str, &suff_arry, read_k);
		//let elapsed_align_initial = start_align_initial.elapsed();
		//println!("Align Initial Time: {} seconds", elapsed_align_initial.as_millis());
		let mut mmp = 0;
		if beg <= end {
			//println!("Interval found: {}",read_k);
			//If an interval is found begin MMP
			let mut beg_p = beg;
			let mut end_p = end;
			
			//finding MMP
			//println!("SA index Match: Beg: {} & End: {}",beg,end);
			let mut done = false; 
			let start_align_mmp = Instant::now();
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
			//let elapsed_align_mmp = start_align_mmp.elapsed();
			//println!("Align MMP Time: {} seconds", elapsed_align_mmp.as_millis());
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
		for (transcript, coverage) in cov_consensus.iter() {
			if *coverage > threshold {
				trans.push((*transcript,*coverage));
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

	//let elapsed_align = start_align.elapsed();
	//println!("Align Time: {} seconds", elapsed_align.as_millis());
	

	trans
}