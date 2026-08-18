#ifndef SP_RUNTIME_H
#define SP_RUNTIME_H

#define CLIP_IF_HIGHER(x, upperBound) if (##x > ##upperBound) ##x = ##upperBound
#define CLIP_IF_LOWER(x, upperBound) if (##x < ##upperBound) ##x = ##upperBound


#endif//!SP_RUNTIME_H
