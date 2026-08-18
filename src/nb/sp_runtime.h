#ifndef SP_RUNTIME_H
#define SP_RUNTIME_H

#define CLIP_IF_HIGHER(x, upperBound) if (##x > ##upperBound) ##x = ##upperBound;
#define CLIP_IF_LOWER(x, upperBound) if (##x < ##upperBound) ##x = ##upperBound;

#define CLIP_IF_LOWER_IDX(x, upperBound, idxTarget, idxData) if (##x < ##upperBound) {\
    ##x = ##upperBound; \
    ##idxTarget = ##idxData; \
}

#define CLIP_IF_HIGHER_IDX(x, upperBound, idxTarget, idxData) if (##x > ##upperBound) {\
    ##x = ##upperBound; \
    ##idxTarget = ##idxData; \
}

#define CLIP_IF_GTQ_IDX(x, upperBound, idxTarget, idxData) if (##x >= ##upperBound) {\
    ##x = ##upperBound; \
    ##idxTarget = ##idxData; \
}

#define SWAP_IF_LOWER(x, upperBound) if (##x < ##upperBound) ##upperBound = ##x;

#define TRY_ALLOC(dest, type) if ( ( ##dest = ( ##type * ) malloc( sizeof( ##type ) ) ) == NULL ) {\
      return-1;\
   }

#endif//!SP_RUNTIME_H
