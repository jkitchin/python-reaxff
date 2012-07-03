import hashlib, pickle
def reax_hash(atoms):
    '''returns the first 15 characters of a hash of the pickled string
    representing an atoms object.

    used to generate unique keys for an atoms object. There is some
    concern that there could be collisions in just the first 15
    characters of the hash, but 15 characters is the limit of what the
    Reax code can handle.

    It is not sufficient to simply make 15 random characters because
    you need to ensure the same 15 characters every time you see a
    particular atoms object.

    This method is fragile to changes in the atoms object and pickling
    method. Any changes to those, e.g. additional information stored
    in the atoms object, or changes in the pickling algorithm, will
    result in changes to the hash. We could replace this with a
    hard-coded string representation of an atoms object and take a
    hash of that.
    '''
    return hashlib.sha224(pickle.dumps(atoms)).hexdigest()[:15]

if __name__ == '__main__':
    print reax_hash('test')
    print reax_hash('testt')
