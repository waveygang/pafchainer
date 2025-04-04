use std::error::Error;

use log::debug;

// CIGAR operation
#[derive(Debug, Clone, Copy)]
pub struct CigarOp(pub char, pub usize);

// Parse CIGAR string into operations
pub fn cigar_to_ops(cigar: &str) -> Result<Vec<CigarOp>, Box<dyn Error>> {
    let mut ops = Vec::new();
    let mut count = 0;

    for c in cigar.chars() {
        if c.is_digit(10) {
            count = count * 10 + c.to_digit(10).unwrap() as usize;
        } else {
            ops.push(CigarOp(c, count));
            count = 0;
        }
    }

    if count > 0 {
        return Err("CIGAR string ended with a number".into());
    }

    Ok(ops)
}

// Convert WFA alignment CIGAR to standard operations
pub fn cigar_u8_to_ops(cigar: &[u8]) -> Vec<CigarOp> {
    if cigar.is_empty() {
        return Vec::new();
    }

    let mut ops = Vec::new();
    let mut count = 0;
    let mut last_op = cigar[0];

    for &op in cigar {
        if op == last_op {
            count += 1;
        } else {
            // Convert the operation
            let cigar_op = match last_op {
                b'M' => '=', // Match
                b'X' => 'X', // Mismatch
                b'D' => 'D', // Deletion
                b'I' => 'I', // Insertion
                _ => panic!("Unknown CIGAR operation: {}", last_op),
            };
            ops.push(CigarOp(cigar_op, count));
            count = 1;
            last_op = op;
        }
    }

    // Add the last operation
    let cigar_op = match last_op {
        b'M' => '=',
        b'X' => 'X',
        b'D' => 'D',
        b'I' => 'I',
        _ => panic!("Unknown CIGAR operation: {}", last_op),
    };
    ops.push(CigarOp(cigar_op, count));
    ops
}

// Convert CIGAR operations to string
pub fn ops_to_cigar(ops: &[CigarOp]) -> String {
    ops.iter()
        .map(|&CigarOp(op, count)| format!("{}{}", count, op))
        .collect()
}

// Erode CIGAR from the end with separate query and target erosion limits
pub fn erode_cigar_end(
    cigar_ops: &[CigarOp],
    query_bp_limit: usize,
    target_bp_limit: usize,
) -> Result<(Vec<CigarOp>, usize, usize), Box<dyn Error>> {
    let mut new_ops_reversed = Vec::new();
    let mut remaining_query_bp = query_bp_limit;
    let mut remaining_target_bp = target_bp_limit;
    let mut query_bases_removed = 0;
    let mut target_bases_removed = 0;

    debug!(
        "\tEroding from end of CIGAR - Query limit: {} bp, Target limit: {} bp",
        query_bp_limit, target_bp_limit
    );

    // Copy operations from beginning until we need to erode
    let mut i = 0;
    while i < cigar_ops.len() && (remaining_query_bp > 0 || remaining_target_bp > 0) {
        // Process operations from the end
        let idx = cigar_ops.len() - 1 - i;
        let CigarOp(op, count) = cigar_ops[idx];

        match op {
            '=' | 'X' | 'M' => {
                // Advances both sequences
                if remaining_query_bp > 0 || remaining_target_bp > 0 {
                    let query_erosion = std::cmp::min(count, remaining_query_bp);
                    let target_erosion = std::cmp::min(count, remaining_target_bp);

                    // Determine actual erosion
                    let actual_erosion = std::cmp::max(query_erosion, target_erosion);

                    if actual_erosion < count {
                        // Partial erosion
                        new_ops_reversed.push(CigarOp(op, count - actual_erosion));
                    }

                    query_bases_removed += actual_erosion;
                    target_bases_removed += actual_erosion;
                    remaining_query_bp = remaining_query_bp.saturating_sub(actual_erosion);
                    remaining_target_bp = remaining_target_bp.saturating_sub(actual_erosion);
                } else {
                    // No more erosion needed
                    new_ops_reversed.push(CigarOp(op, count));
                }
            }
            'I' => {
                // Insertion only advances query
                if remaining_query_bp > 0 {
                    let actual_erosion = std::cmp::min(count, remaining_query_bp);

                    if actual_erosion < count {
                        // Partial erosion
                        new_ops_reversed.push(CigarOp(op, count - actual_erosion));
                    }

                    query_bases_removed += actual_erosion;
                    remaining_query_bp = remaining_query_bp.saturating_sub(actual_erosion);
                } else {
                    // No more query erosion needed
                    new_ops_reversed.push(CigarOp(op, count));
                }
            }
            'D' => {
                // Deletion only advances target
                if remaining_target_bp > 0 {
                    let actual_erosion = std::cmp::min(count, remaining_target_bp);

                    if actual_erosion < count {
                        // Partial erosion
                        new_ops_reversed.push(CigarOp(op, count - actual_erosion));
                    }

                    target_bases_removed += actual_erosion;
                    remaining_target_bp = remaining_target_bp.saturating_sub(actual_erosion);
                } else {
                    // No more target erosion needed
                    new_ops_reversed.push(CigarOp(op, count));
                }
            }
            _ => new_ops_reversed.push(CigarOp(op, count)), // Keep other operations
        }

        i += 1;
    }

    // Reverse the eroded part to restore the original order.
    new_ops_reversed.reverse();

    // Combine the unprocessed (beginning) portion with the processed (eroded) part.
    let mut new_ops = Vec::with_capacity(cigar_ops.len() - i + new_ops_reversed.len());
    new_ops.extend_from_slice(&cigar_ops[..(cigar_ops.len() - i)]);
    new_ops.extend(new_ops_reversed);

    Ok((new_ops, query_bases_removed, target_bases_removed))
}

// Erode CIGAR from the start with separate query and target erosion limits
pub fn erode_cigar_start(
    cigar_ops: &[CigarOp],
    query_bp_limit: usize,
    target_bp_limit: usize,
) -> Result<(Vec<CigarOp>, usize, usize), Box<dyn Error>> {
    let mut new_ops = Vec::new();
    let mut remaining_query_bp = query_bp_limit;
    let mut remaining_target_bp = target_bp_limit;
    let mut query_bases_removed = 0;
    let mut target_bases_removed = 0;

    debug!(
        "\tEroding from start of CIGAR - Query limit: {} bp, Target limit: {} bp",
        query_bp_limit, target_bp_limit
    );

    // Skip or partially include first operations until eroded enough
    let mut i = 0;
    while i < cigar_ops.len() && (remaining_query_bp > 0 || remaining_target_bp > 0) {
        let CigarOp(op, count) = cigar_ops[i];

        match op {
            '=' | 'X' | 'M' => {
                // Advances both sequences
                if remaining_query_bp > 0 || remaining_target_bp > 0 {
                    let query_erosion = std::cmp::min(count, remaining_query_bp);
                    let target_erosion = std::cmp::min(count, remaining_target_bp);

                    // Determine actual erosion
                    let actual_erosion = std::cmp::max(query_erosion, target_erosion);

                    if actual_erosion < count {
                        // Partial erosion
                        new_ops.push(CigarOp(op, count - actual_erosion));
                    }

                    query_bases_removed += actual_erosion;
                    target_bases_removed += actual_erosion;
                    remaining_query_bp = remaining_query_bp.saturating_sub(actual_erosion);
                    remaining_target_bp = remaining_target_bp.saturating_sub(actual_erosion);
                } else {
                    // No more erosion needed
                    new_ops.push(CigarOp(op, count));
                }
            }
            'I' => {
                // Insertion only advances query
                if remaining_query_bp > 0 {
                    let actual_erosion = std::cmp::min(count, remaining_query_bp);

                    if actual_erosion < count {
                        // Partial erosion
                        new_ops.push(CigarOp(op, count - actual_erosion));
                    }

                    query_bases_removed += actual_erosion;
                    remaining_query_bp = remaining_query_bp.saturating_sub(actual_erosion);
                } else {
                    // No more query erosion needed
                    new_ops.push(CigarOp(op, count));
                }
            }
            'D' => {
                // Deletion only advances target
                if remaining_target_bp > 0 {
                    let actual_erosion = std::cmp::min(count, remaining_target_bp);

                    if actual_erosion < count {
                        // Partial erosion
                        new_ops.push(CigarOp(op, count - actual_erosion));
                    }

                    target_bases_removed += actual_erosion;
                    remaining_target_bp = remaining_target_bp.saturating_sub(actual_erosion);
                } else {
                    // No more target erosion needed
                    new_ops.push(CigarOp(op, count));
                }
            }
            _ => new_ops.push(CigarOp(op, count)), // Keep other operations
        }

        i += 1;
    }

    // Add all remaining operations
    while i < cigar_ops.len() {
        new_ops.push(cigar_ops[i]);
        i += 1;
    }

    Ok((new_ops, query_bases_removed, target_bases_removed))
}

// Merge multiple CIGAR operations vectors
pub fn merge_cigar_ops(cigar_ops_list: &[&[CigarOp]]) -> Result<Vec<CigarOp>, Box<dyn Error>> {
    // Estimate capacity to avoid reallocations
    let total_capacity: usize = cigar_ops_list.iter().map(|ops| ops.len()).sum();
    let mut optimized_ops = Vec::with_capacity(total_capacity);

    // Track current operation type
    let mut current_op: Option<CigarOp> = None;

    // Process all operations from all input vectors sequentially
    for &cigar_ops in cigar_ops_list {
        for &op in cigar_ops {
            if let Some(CigarOp(prev_op, prev_count)) = current_op {
                if prev_op == op.0 {
                    // Combine with previous operation
                    current_op = Some(CigarOp(prev_op, prev_count + op.1));
                } else {
                    // Different operation, push previous and update current
                    optimized_ops.push(CigarOp(prev_op, prev_count));
                    current_op = Some(op);
                }
            } else {
                // First operation
                current_op = Some(op);
            }
        }
    }

    // Don't forget to add the last operation
    if let Some(op) = current_op {
        optimized_ops.push(op);
    }

    Ok(optimized_ops)
}

// Calculate CIGAR statistics
pub fn calculate_cigar_stats(
    cigar_ops: &[CigarOp],
) -> Result<(u64, u64, u64, u64, u64, u64, u64), Box<dyn Error>> {
    let (matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, block_len) =
        cigar_ops.iter().fold(
            (0, 0, 0, 0, 0, 0, 0),
            |(m, mm, i, i_bp, d, d_bp, bl), &CigarOp(op, len)| {
                let len = len as u64;
                match op {
                    'M' => (m + len, mm, i, i_bp, d, d_bp, bl + len),
                    '=' => (m + len, mm, i, i_bp, d, d_bp, bl + len),
                    'X' => (m, mm + len, i, i_bp, d, d_bp, bl + len),
                    'I' => (m, mm, i + 1, i_bp + len, d, d_bp, bl + len),
                    'D' => (m, mm, i, i_bp, d + 1, d_bp + len, bl + len),
                    _ => (m, mm, i, i_bp, d, d_bp, bl),
                }
            },
        );

    Ok((
        matches,
        mismatches,
        insertions,
        inserted_bp,
        deletions,
        deleted_bp,
        block_len,
    ))
}
