/* eslint-disable no-underscore-dangle */
/* eslint-disable react/no-array-index-key */

import React from 'react'
import PropTypes from 'prop-types'

import { Confirm, Form } from 'semantic-ui-react'
import isEqual from 'lodash/isEqual'

import { HttpRequestHelper } from '../../utils/httpRequestHelper'
import RequestStatus from '../form/RequestStatus'
import ButtonPanel from './ButtonPanel'
import MessagesPanel from './MessagesPanel'

/**
 * Form wrapper that provides Submit and Cancel functionality.
 */
class FormWrapper extends React.Component
{
  static propTypes = {
    cancelButtonText: PropTypes.string,
    submitButtonText: PropTypes.string,
    getFormDataJson: PropTypes.func, // required if either formSubmitUrl or performValidation is provided
    formSubmitUrl: PropTypes.string,
    performClientSideValidation: PropTypes.func,
    handleSave: PropTypes.func,
    handleClose: PropTypes.func,
    confirmCloseIfNotSaved: PropTypes.bool.isRequired,
    size: PropTypes.string, // form size (see https://react.semantic-ui.com/collections/form#form-example-size)
    children: PropTypes.node,
  }

  constructor(props) {
    super(props)

    this.preventSubmit = null

    this.state = {
      saveStatus: RequestStatus.NONE, // one of NONE, IN_PROGRESS, SUCCEEDED, ERROR
      saveErrorMessage: null,
      isConfirmCloseVisible: false, // show confirm dialog to ask user if they really want to close without saving
      errors: [],
      warnings: [],
      info: [],
    }

    if (this.props.confirmCloseIfNotSaved) {
      this.originalFormData = JSON.parse(JSON.stringify(this.props.getFormDataJson())) //make a deep copy
    }
  }

  doClientSideValidation() {
    if (!this.props.performClientSideValidation) {
      return
    }

    let validationResult = this.props.performClientSideValidation(this.props.getFormDataJson())
    if (validationResult === null) {
      validationResult = {}
    }

    this.setState({
      errors: validationResult.errors,
      warnings: validationResult.warnings,
      info: validationResult.info,
    })

    if (validationResult.preventSubmit || (validationResult.errors && validationResult.errors.length > 0)) {
      this.preventSubmit = 'clientSideValidation'
    } else {
      this.preventSubmit = null
    }
  }

  /*
  componentWillReceiveProps() {
    this.doClientSideValidation()
  }
  */

  componentWillMount() {
    this.setState({
      saveStatus: RequestStatus.NONE,
      saveErrorMessage: null,
    })
  }

  hasFormBeenModified = () => {
    console.log('this.originalFormData', this.originalFormData)
    console.log('this.props.getFormDataJson()', this.props.getFormDataJson())
    return !isEqual(this.originalFormData, this.props.getFormDataJson())
  }

  doSave = (e) => {
    e.preventDefault()

    //do client-side validation
    this.doClientSideValidation()

    if (this.preventSubmit === 'clientSideValidation') {
      return // don't submit the form if there are error
    }

    this.setState({
      saveStatus: RequestStatus.IN_PROGRESS,
      saveErrorMessage: null,
    })

    if (this.props.formSubmitUrl) {
      const httpRequestHelper = new HttpRequestHelper(
        this.props.formSubmitUrl,
        this.handleFormSubmitSucceeded,
        this.handleFormSubmitError,
      )

      const jsonContents = {
        form: this.props.getFormDataJson(),
      }

      //console.log('Posting: ', jsonContents)
      httpRequestHelper.post(jsonContents)
    } else {
      this.doClose(false)
    }
  }

  doClose = (confirmCloseIfNecessary) => {
    if (confirmCloseIfNecessary && this.props.confirmCloseIfNotSaved && this.hasFormBeenModified()) {
      //first double check that user wants to close
      this.setState({
        isConfirmCloseVisible: true,
      })
    } else if (this.props.handleClose) {
      this.props.handleClose()
    }
  }

  handleFormSubmitSucceeded = (responseJson) => {

    console.log('got response: ', responseJson)
    //allow server-side validation
    if (responseJson.errors) {
      this.setState({
        saveStatus: RequestStatus.NONE,
        errors: responseJson.errors,
      })
      this.preventSubmit = 'serverSideValidation'
      return //cancel submit
    }

    this.setState({
      saveStatus: RequestStatus.NONE,
      saveErrorMessage: null,
    })

    if (this.props.handleSave) {
      this.props.handleSave(responseJson)
    }
    this.doClose(false)
  }

  handleFormSubmitError = (exception) => {
    console.error('Exception in HttpRequestHelper:', exception)
    //if saveStatus === RequestStatus.NONE, the component has been reset
    if (this.state.saveStatus !== RequestStatus.NONE) {
      this.setState({
        saveStatus: RequestStatus.ERROR,
        saveErrorMessage: exception.message.toString(),
      })
    }
  }

  renderForm() {
    const children = React.Children.map(
      this.props.children,
      child => (
        child.props && child.props.name && this.state.errors && this.state.errors[child.props.name] !== undefined ?
          React.cloneElement(child, { error: true }) : child
      ),
    )

    return (
      <Form onSubmit={this.doSave} style={{ textAlign: 'left' }} size={this.props.size}>
        {children}
      </Form>
    )
  }

  renderMessageBoxes() {
    return <MessagesPanel errors={this.state.errors} warnings={this.state.warnings} info={this.state.info} />
  }

  renderButtonPanel() {
    return <ButtonPanel
      cancelButtonText={this.props.cancelButtonText}
      submitButtonText={this.props.submitButtonText}
      handleClose={(e) => { e.preventDefault(); this.doClose(true) }}
      handleSave={this.doSave}
      saveStatus={this.state.saveStatus}
      saveErrorMessage={this.state.saveErrorMessage}
    />
  }

  renderConfirmCloseDialog() {
    return <Confirm
      content="The form contains unsaved changes. Are you sure you want to close it?"
      open={this.state.isConfirmCloseVisible}
      onCancel={() => this.setState({ isConfirmCloseVisible: false })}
      onConfirm={() => this.doClose(false)}
    />
  }

  render() {
    return (
      <div>
        {this.renderForm()}
        {this.renderMessageBoxes()}
        {this.renderButtonPanel()}
        {this.renderConfirmCloseDialog()}
      </div>
    )
  }

}


export default FormWrapper
