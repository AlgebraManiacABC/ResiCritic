class Peptide:
    def __init__(self,
                 seq: str = "",
                 mod_list: list[str] = None,
                 counts: dict[str,int] = None):

        self.sequence = seq
        self.mod_list = mod_list
        self.counts = counts

    def __str__(self):
        # Output as follows:
        # ##:ABCDE*FG<#;#;#>
        # ## == Starting index (1-indexed)
        # ABCDE*FG == Sequence, with modifications included
        # <#;#;#> == Counts per column
        local_ids = [int(resi[1:]) for resi in self.mod_list]
        string = self.sequence
        for i in sorted(local_ids,reverse=True):
            string = string[:i] + '*' + string[i:]
        string += "<" + ''.join([f"{count};" for count in self.counts.values()[:-1] + ">"])
        return string

    def to_str_by_fasta(self, fasta_sequence, residue_id: int=0):
        zero_index = fasta_sequence.find(self.sequence)
        if zero_index < 0:
            return ""

        local_ids = [int(resi[1:]) for resi in self.mod_list]
        string = self.sequence
        for i in range(len(self.sequence), 0, -1):
            if i in local_ids:
                string = string[:i] + '*' + string[i:]
            if i+zero_index == residue_id:
                string = string[:i-1] + '(' + string[i:i+2] + ')' + string[i+2:]
        string = f"{zero_index+1}:" + string
        string += "<" + ''.join([f"{count};" for count in self.counts.values()[:-1] + ">"])
        return string
