import { combineReducers } from 'redux'
import { SubmissionError } from 'redux-form'

import { loadingReducer, createSingleObjectReducer, createObjectsByIdReducer, createSingleValueReducer } from 'redux/utils/reducerFactories'
import { REQUEST_PROJECTS, RECEIVE_DATA } from 'redux/rootReducer'
import { HttpRequestHelper } from 'shared/utils/httpRequestHelper'
import { getProject, getProjectFamilies } from './selectors'
import { SHOW_ALL, SORT_BY_FAMILY_NAME, SORT_BY_FAMILY_GUID } from './constants'

// action creators and reducers in one file as suggested by https://github.com/erikras/ducks-modular-redux

const UPDATE_FAMILY_TABLE_STATE = 'UPDATE_FAMILY_TABLE_STATE'
const UPDATE_SAVED_VARIANT_TABLE_STATE = 'UPDATE_VARIANT_STATE'
const UPDATE_CURRENT_PROJECT = 'UPDATE_CURRENT_PROJECT'
const REQUEST_PROJECT_DETAILS = 'REQUEST_PROJECT_DETAILS'
const RECEIVE_PROJECT_DETAILS = 'RECEIVE_PROJECT_DETAILS'
const REQUEST_SAVED_VARIANTS = 'REQUEST_SAVED_VARIANTS'
const RECEIVE_SAVED_VARIANTS = 'RECEIVE_SAVED_VARIANTS'
const RECEIVE_SAVED_VARIANT_FAMILIES = 'RECEIVE_SAVED_VARIANT_FAMILIES'

// Data actions

export const loadProject = (projectGuid) => {
  return (dispatch, getState) => {
    dispatch({ type: UPDATE_CURRENT_PROJECT, newValue: projectGuid })

    const currentProject = getState().projectsByGuid[projectGuid]
    if (!currentProject || !currentProject.detailsLoaded) {
      dispatch({ type: REQUEST_PROJECT_DETAILS })
      if (!currentProject) {
        dispatch({ type: REQUEST_PROJECTS })
      }
      new HttpRequestHelper(`/api/project/${projectGuid}/details`,
        (responseJson) => {
          dispatch({ type: RECEIVE_PROJECT_DETAILS })
          dispatch({ type: RECEIVE_DATA, updatesById: { projectsByGuid: { [projectGuid]: responseJson.project }, ...responseJson } })
        },
        (e) => {
          dispatch({ type: RECEIVE_DATA, error: e.message, updatesById: {} })
        },
      ).get()
    }
  }
}

export const loadProjectVariants = (familyGuid, variantGuid) => {
  return (dispatch, getState) => {
    const state = getState()
    const project = getProject(state)

    let url = `/api/project/${project.projectGuid}/saved_variants`

    // Do not load if already loaded
    let expectedFamilyGuids
    if (variantGuid) {
      if (state.projectSavedVariants[variantGuid]) {
        return
      }
      url = `${url}/${variantGuid}`
    } else {
      expectedFamilyGuids = familyGuid ? [familyGuid] : getProjectFamilies(state).map(family => family.familyGuid)
      if (expectedFamilyGuids.length > 0 && expectedFamilyGuids.every(family => state.projectSavedVariantFamilies[family])) {
        return
      }
    }

    dispatch({ type: REQUEST_SAVED_VARIANTS })
    new HttpRequestHelper(url,
      (responseJson) => {
        dispatch({ type: RECEIVE_SAVED_VARIANTS, updatesById: responseJson.savedVariants })
        if (expectedFamilyGuids) {
          dispatch({
            type: RECEIVE_SAVED_VARIANT_FAMILIES,
            updates: expectedFamilyGuids.reduce((acc, family) => ({ ...acc, [family]: true }), {}),
          })
        }
      },
      (e) => {
        dispatch({ type: RECEIVE_SAVED_VARIANTS, error: e.message, updatesById: {} })
      },
    ).get(familyGuid ? { family: familyGuid } : {})
  }
}

export const unloadProject = () => {
  return (dispatch, getState) => {
    const state = getState()
    const deleteVariants = Object.keys(state.projectSavedVariants).reduce((acc, o) => ({ ...acc, [o]: null }))
    const deleteVariantFamilies = Object.keys(state.projectSavedVariantFamilies).reduce((acc, o) => ({ ...acc, [o]: false }))
    dispatch({ type: UPDATE_CURRENT_PROJECT, newValue: null })
    dispatch({ type: REQUEST_SAVED_VARIANTS, updatesById: deleteVariants })
    dispatch({ type: RECEIVE_SAVED_VARIANT_FAMILIES, updates: deleteVariantFamilies })
  }
}


export const updateFamilies = (values) => {
  return (dispatch, getState) => {
    const action = values.delete ? 'delete' : 'edit'
    return new HttpRequestHelper(`/api/project/${getState().currentProjectGuid}/${action}_families`,
      (responseJson) => {
        dispatch({ type: RECEIVE_DATA, updatesById: responseJson })
      },
      (e) => { throw new SubmissionError({ _error: [e.message] }) },
    ).post(values)
  }
}

export const updateIndividuals = (values) => {
  return (dispatch, getState) => {
    let action = 'edit_individuals'
    if (values.uploadedFileId) {
      action = `save_individuals_table/${values.uploadedFileId}`
    } else if (values.delete) {
      action = 'delete_individuals'
    }

    return new HttpRequestHelper(`/api/project/${getState().currentProjectGuid}/${action}`,
      (responseJson) => {
        dispatch({ type: RECEIVE_DATA, updatesById: responseJson })
      },
      (e) => {
        if (e.body && e.body.errors) {
          throw new SubmissionError({ _error: e.body.errors })
          // e.body.warnings.forEach((err) => { throw new SubmissionError({ _warning: err }) })
        } else {
          throw new SubmissionError({ _error: [e.message] })
        }
      },
    ).post(values)
  }
}

// Family table actions

export const setCurrentPage = currentPage => ({ type: UPDATE_FAMILY_TABLE_STATE, updates: { currentPage } })
export const setRecordsPerPage = recordsPerPage => ({ type: UPDATE_FAMILY_TABLE_STATE, updates: { recordsPerPage } })

export const updateFamiliesFilter = familiesFilter => ({ type: UPDATE_FAMILY_TABLE_STATE, updates: { familiesFilter } })
export const updateFamiliesSortOrder = familiesSortOrder => ({ type: UPDATE_FAMILY_TABLE_STATE, updates: { familiesSortOrder } })
export const updateFamiliesSortDirection = familiesSortDirection => ({ type: UPDATE_FAMILY_TABLE_STATE, updates: { familiesSortDirection } })
export const updateShowDetails = showDetails => ({ type: UPDATE_FAMILY_TABLE_STATE, updates: { showDetails } })

export const updateSavedVariantTable = updates => ({ type: UPDATE_SAVED_VARIANT_TABLE_STATE, updates })

// reducers

export const reducers = {
  currentProjectGuid: createSingleValueReducer(UPDATE_CURRENT_PROJECT, null),
  projectDetailsLoading: loadingReducer(REQUEST_PROJECT_DETAILS, RECEIVE_PROJECT_DETAILS),
  projectSavedVariants: createObjectsByIdReducer(RECEIVE_SAVED_VARIANTS),
  projectSavedVariantsLoading: loadingReducer(REQUEST_SAVED_VARIANTS, RECEIVE_SAVED_VARIANTS),
  projectSavedVariantFamilies: createSingleObjectReducer(RECEIVE_SAVED_VARIANT_FAMILIES),
  familyTableState: createSingleObjectReducer(UPDATE_FAMILY_TABLE_STATE, {
    currentPage: 1,
    recordsPerPage: 100,
    familiesFilter: SHOW_ALL,
    familiesSortOrder: SORT_BY_FAMILY_NAME,
    familiesSortDirection: 1,
    showDetails: true,
  }, false),
  savedVariantTableState: createSingleObjectReducer(UPDATE_SAVED_VARIANT_TABLE_STATE, {
    hideExcluded: false,
    categoryFilter: SHOW_ALL,
    sortOrder: SORT_BY_FAMILY_GUID,
    currentPage: 1,
    recordsPerPage: 25,
  }, false),
}

const rootReducer = combineReducers(reducers)

export default rootReducer