"""Version 3: v2 preparation with independent mutable model metadata.

The frozen v2 algorithm, policy types and report schema are reused unchanged.
COBRA's fast Model.copy() retains references to some nested metadata and to
compartment descriptions. Running preparation on a deep copy detaches that
state, including groups and solver objects, from the caller before preparation.
"""
from __future__ import annotations

from copy import deepcopy
from typing import Iterable

import cobra

from .media import Medium
from .medium_preparation_v2 import CompartmentPolicy, PreparedMedium
from .medium_preparation_v2 import prepare_medium as _prepare_medium_v2


def prepare_medium(model: cobra.Model, medium: Medium, *, policy: CompartmentPolicy,
                   completion_media: Iterable[Medium] = (), missing_policy: str = 'error') -> PreparedMedium:
    """Return v2's prepared formulation and report with a fully detached model.

    All exchange validation, explicit compartment handling, completion rules and
    missing-component behavior are inherited from v2. Preparation still precedes
    gene deletions and does not optimize a model. Nested metadata, compartment
    descriptions, groups and solver state can be edited independently in either
    the returned model or the caller's model. The caller's active contexts remain
    intact; the returned model has its own empty context stack.

    The additional deep copy must succeed: a copy failure propagates without
    returning a partially isolated model or changing the caller's input.
    """
    # COBRA clears contexts during deepcopy. Copy before invoking v2 so that
    # fast-copy bookkeeping cannot register undo callbacks on a caller context.
    # V2's copy also rebuilds group membership and group-to-model ownership,
    # which COBRA's serialization-based deepcopy alone does not restore.
    isolated = deepcopy(model)
    return _prepare_medium_v2(isolated, medium, policy=policy,
                              completion_media=completion_media, missing_policy=missing_policy)
