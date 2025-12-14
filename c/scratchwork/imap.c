struct imap_s
{
	int regions;
	int *indices;
	int offset;
	int length;
};

int imap_init(struct imap_s *imap, int regions, int length)
{
	int alloc_status = 0;
	
	imap->indices = calloc(2 * regions + length, sizeof(int));
	
	if (imap->indices) {
		imap->length = length;
		imap->offset = 2 * regions;
	}
	
	return imap->indices != NULL;
}

void imap_set_region_bounds(struct imap_s *imap, int region, int start, int end)
{
	imap->indices[2 * region] = start;
	imap->indices[2 * region] = end;
}

void imap_set_region_indices(struct imap_s *imap, int region, int istart, int iend, int jstart, int jend)
{
	int *indices = imap->indices + 2*region;
	int i, j, position = 0;
	
	for (j = jstart; j <= jend; j++) {
		jind = j * M;
		
		for (i = istart; i <= iend; i++) {
			indices[position++] = jind + i;
		}
	}
}